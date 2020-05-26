/**
 * @file    src/pem_bwt_src/compute_bwt.hpp
 * @section LICENCE
 *
 * This file is part of pEM-BWT v0.1.0
 * Copyright (C) 2016-2020
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __SRC_PEM_BWT_SRC_COMPUTE_BWT_HPP_INCLUDED
#define __SRC_PEM_BWT_SRC_COMPUTE_BWT_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <string>
#include <algorithm>
#include <omp.h>
#include <unistd.h>

#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_multi_stream_reader.hpp"
#include "io/async_multi_stream_writer.hpp"
#include "utils.hpp"


namespace pem_bwt_private {

template<typename char_type,
  typename text_offset_type>
void compute_bwt(
    std::string text_filename,
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t ram_use) {

  // Initialize utils and PRNG.
  utils::initialize_stats();
  srand(time(0) + getpid());

  // Empty page cache.
  utils::empty_page_cache(text_filename);
  utils::empty_page_cache(sa_filename);

  // Start the timer.
  long double global_start = utils::wclock();

  // Compute basic parameters.
  std::uint64_t text_length =
    utils::file_size(text_filename) /
    sizeof(char_type);
  std::uint64_t total_io_volume = 0;

  // Set optimal buffer sizes.
  static const std::uint64_t opt_io_buf_ram = ((std::uint64_t)8 << 20);
  static const std::uint64_t opt_local_buf_ram = ((std::uint64_t)8 << 20);
  std::uint64_t max_block_size_ram =
    (std::uint64_t)((long double)ram_use * 0.9L);
  std::uint64_t io_buf_ram = opt_io_buf_ram;
  std::uint64_t local_buf_ram = opt_local_buf_ram;

  // Shrink buffers if necessary, otherwise extend block.
  if (max_block_size_ram + io_buf_ram + local_buf_ram > ram_use) {
    std::uint64_t total_buf_ram = io_buf_ram + local_buf_ram;
    std::uint64_t ram_budget = ram_use - max_block_size_ram;
    long double shrink_factor =
      ((long double)ram_budget) / ((long double)total_buf_ram);
    io_buf_ram = (std::uint64_t)((long double)io_buf_ram * shrink_factor);
    local_buf_ram =
      (std::uint64_t)((long double)local_buf_ram * shrink_factor);
  } else max_block_size_ram = ram_use - io_buf_ram - local_buf_ram;

  // Compute block size.
  std::uint64_t max_block_size =
    utils::disk_block_size<char_type>(max_block_size_ram);
  std::uint64_t n_blocks =
    (text_length + max_block_size - 1) / max_block_size;

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  sa_filename = utils::absolute_path(sa_filename);
  output_filename = utils::absolute_path(output_filename);

  // Print summary of basic parameters.
  fprintf(stderr, "Running pEM-BWT v0.1.0\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "SA filename = %s\n", sa_filename.c_str());
  fprintf(stderr, "Output (BWT) filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n",
      text_length, (1.L * text_length * sizeof(char_type)) / (1 << 20));
  fprintf(stderr, "RAM use = %lu bytes (%.2LfMiB)\n",
      ram_use, (1.L * ram_use) / (1 << 20));
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));
  fprintf(stderr, "Parallel mode = ");

#ifdef _OPENMP
  fprintf(stderr, "ON\n");
  fprintf(stderr, "Number of threads = %d\n", omp_get_max_threads());
#else
  fprintf(stderr, "OFF\n");
#endif

  fprintf(stderr, "\n\n");

  fprintf(stderr, "Compute BWT from text and SA:\n");
  fprintf(stderr, "  Segment size = %lu symbols (%.2LfMiB)\n",
      max_block_size,
      (1.L * max_block_size * sizeof(char_type)) / (1L << 20));
  fprintf(stderr, "  Number of segments = %lu\n", n_blocks);

  // Set the filenames of files storing SA and BWT subsequences.
  std::string *sa_subsequences_filenames = new std::string[n_blocks];
  std::string *bwt_subsequences_filenames = new std::string[n_blocks];
  for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
    sa_subsequences_filenames[block_id] =
      output_filename + ".sa_subseq." +
      utils::intToStr(block_id) + "." +
      utils::random_string_hash();
    bwt_subsequences_filenames[block_id] =
      output_filename + ".bwt_sebseq." +
      utils::intToStr(block_id) + "." +
      utils::random_string_hash();
  }

  // Compute SA subsequences.
  if (n_blocks > 1) {

    // Print initial message and start the timer.
    fprintf(stderr, "  Compute SA subsequences: ");
    long double compute_sa_subseq_start = utils::wclock();

    // Initialize I/O volume.
    std::uint64_t io_volume = 0;

    // Initialize reader of SA.
    typedef async_stream_reader<text_offset_type> sa_reader_type;
    sa_reader_type *sa_reader = new sa_reader_type(sa_filename, io_buf_ram,
        std::max((std::uint64_t)4, io_buf_ram / ((std::uint64_t)2 << 20)));

    // Initialize writer of SA subsequences.
    static const std::uint64_t n_empty_buffers = 4;
    std::uint64_t n_total_bufs = n_blocks + n_empty_buffers;
    std::uint64_t total_buffers_ram = max_block_size_ram;
    std::uint64_t buffer_size =
      std::min((16UL << 20), total_buffers_ram / n_total_bufs);
    typedef async_multi_stream_writer<text_offset_type> sa_subseq_writer_type;
    sa_subseq_writer_type *sa_subseq_writer =
      new sa_subseq_writer_type(n_blocks, buffer_size, n_empty_buffers);
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id)
      sa_subseq_writer->add_file(sa_subsequences_filenames[block_id]);

    // Read SA / write SA subsequences.
    for (std::uint64_t j = 0; j < text_length; ++j) {
      std::uint64_t sa_j = sa_reader->read();
      if (sa_j == 0) continue;
      else --sa_j;
      std::uint64_t block_id = sa_j / max_block_size;
      sa_subseq_writer->write_to_ith_file(block_id, sa_j);
    }

    // Stop I/O threads.
    sa_reader->stop_reading();

    // Update I/O volume.
    io_volume +=
      sa_reader->bytes_read() +
      sa_subseq_writer->bytes_written();
    total_io_volume += io_volume;

    // Clean up.
    delete sa_subseq_writer;
    delete sa_reader;

    // Print summary.
    long double compute_sa_subseq_time =
      utils::wclock() - compute_sa_subseq_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, "
        "total I/O = %.2Lfbytes/symbol\n",
        compute_sa_subseq_time,
        ((1.L * io_volume) / (1L << 20)) / compute_sa_subseq_time,
        (1.L * total_io_volume) / text_length);
  }

  // Compute BWT subsequences.
  {

    // Print the initial message and start the timer.
    fprintf(stderr, "  Compute BWT subsequences: ");
    long double compute_bwt_subseq_start = utils::wclock();

    // Initialize I/O volume.
    std::uint64_t io_volume = 0;

    // Allocate the array holding the block of text.
    char_type *text_block =
      utils::allocate_array<char_type>(max_block_size);

    // Compute local buffer size.
#ifndef _OPENMP
    std::uint64_t local_buf_size = local_buf_ram /
      (1 * sizeof(text_offset_type));
#else
    std::uint64_t local_buf_size = local_buf_ram /
      (1 * sizeof(text_offset_type) +
       1 * sizeof(char_type));
#endif

    local_buf_size =
      std::max(local_buf_size, (std::uint64_t)1);

    // Allocate local buffers.
    text_offset_type *sa_buffer =
      utils::allocate_array<text_offset_type>(local_buf_size);

#ifdef _OPENMP
    char_type *bwt_buffer =
      utils::allocate_array<char_type>(local_buf_size);
#endif

    // Process blocks left to right.
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      std::uint64_t block_beg = block_id * max_block_size;
      std::uint64_t block_end =
        std::min(block_beg + max_block_size, text_length);
      std::uint64_t block_size = block_end - block_beg;
      
      // Read block of text into RAM.
      utils::read_at_offset(text_block,
          block_beg * sizeof(char_type),
          block_size, text_filename);
      io_volume += block_size * sizeof(char_type);

      // Compute BWT subsequence and write to file.
      {

        // Initialize SA subsequence reader.
        typedef async_stream_reader<text_offset_type>
          sa_subseq_reader_type;
        std::string sa_subseq_filename =
          (n_blocks == 1) ? sa_filename :
          sa_subsequences_filenames[block_id];
        sa_subseq_reader_type *sa_subseq_reader =
          new sa_subseq_reader_type(sa_subseq_filename,
              io_buf_ram / 2, std::max((std::uint64_t)4,
                (io_buf_ram / 2) / ((std::uint64_t)2 << 20)));

        // Initialize BWT subsequence writer.
        typedef async_stream_writer<char_type>
          bwt_subseq_writer_type;
        std::string bwt_subseq_filename =
          (n_blocks == 1) ? output_filename :
          bwt_subsequences_filenames[block_id];
        bwt_subseq_writer_type *bwt_subseq_writer =
          new bwt_subseq_writer_type(bwt_subseq_filename,
            io_buf_ram / 2, std::max((std::uint64_t)4,
              (io_buf_ram / 2) / ((std::uint64_t)2 << 20)));

        // Compute BWT subsequence.
        std::uint64_t subseq_size =
          utils::file_size(sa_subseq_filename) / sizeof(text_offset_type);
        if (n_blocks == 1) {
          std::uint64_t items_processed = 0;
          while (items_processed < subseq_size) {
            std::uint64_t filled =
              std::min(local_buf_size, subseq_size - items_processed);
            sa_subseq_reader->read(sa_buffer, filled);

#ifdef _OPENMP
            #pragma omp parallel for
            for (std::uint64_t j = 0; j < filled; ++j) {
              std::uint64_t sa_val = sa_buffer[j];
              if (sa_val == 0) bwt_buffer[j] = 0;
              else bwt_buffer[j] = text_block[sa_val - 1];
            }

            bwt_subseq_writer->write(bwt_buffer, filled);
#else
            for (std::uint64_t j = 0; j < filled; ++j) {
              std::uint64_t sa_val = sa_buffer[j];
              char_type bwt_val = 0;
              if (sa_val == 0) bwt_val = 0;
              else bwt_val = text_block[sa_val - 1];
              bwt_subseq_writer->write(bwt_val);
            }
#endif

            items_processed += filled;
          }
        } else {
          std::uint64_t items_processed = 0;
          while (items_processed < subseq_size) {
            std::uint64_t filled =
              std::min(local_buf_size, subseq_size - items_processed);
            sa_subseq_reader->read(sa_buffer, filled);

#ifdef _OPENMP
            #pragma omp parallel for
            for (std::uint64_t j = 0; j < filled; ++j) {
              std::uint64_t sa_val = sa_buffer[j];
              bwt_buffer[j] = text_block[sa_val - block_beg];
            }

            bwt_subseq_writer->write(bwt_buffer, filled);
#else
            for (std::uint64_t j = 0; j < filled; ++j) {
              std::uint64_t sa_val = sa_buffer[j];
              char_type ch = text_block[sa_val - block_beg];
              bwt_subseq_writer->write(ch);
            }
#endif

            items_processed += filled;
          }
        }

        // Stop I/O threads.
        sa_subseq_reader->stop_reading();

        // Update I/O volume.
        io_volume +=
          sa_subseq_reader->bytes_read() +
          bwt_subseq_writer->bytes_written();

        // Clean up.
        delete bwt_subseq_writer;
        delete sa_subseq_reader;
      }

      if (n_blocks > 1)
        utils::file_delete(sa_subsequences_filenames[block_id]);
    }

    // Update I/O volume.
    total_io_volume += io_volume;

    // Clean up.
#ifdef _OPENMP
    utils::deallocate(bwt_buffer);
#endif

    utils::deallocate(sa_buffer);
    utils::deallocate(text_block);

    // Print summary.
    long double compute_bwt_subseq_time =
      utils::wclock() - compute_bwt_subseq_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, "
        "total I/O = %.2Lfbytes/symbol\n",
        compute_bwt_subseq_time,
        ((1.L * io_volume) / (1L << 20)) / compute_bwt_subseq_time,
        (1.L * total_io_volume) / text_length);
  }

  // Merge BWT subsequences.
  if (n_blocks > 1) {

    // Print initial message and start the timer.
    fprintf(stderr, "  Merge BWT subsequences: ");
    long double merge_bwt_subseq_start = utils::wclock();

    // Initialize I/O volume.
    std::uint64_t io_volume = 0;

    // Initialize the reader of BWT subsequences.
    std::uint64_t total_buffers_ram = max_block_size_ram;
    std::uint64_t buffer_size =
      std::max((std::uint64_t)1, total_buffers_ram / n_blocks);
    typedef async_multi_stream_reader<char_type> bwt_subseq_reader_type;
    bwt_subseq_reader_type *bwt_subseq_reader =
      new bwt_subseq_reader_type(n_blocks, buffer_size);
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id)
      bwt_subseq_reader->add_file(bwt_subsequences_filenames[block_id]);

    // Initialize the writer of the final BWT.
    typedef async_stream_writer<char_type> bwt_writer_type;
    bwt_writer_type *bwt_writer = new bwt_writer_type(output_filename,
        io_buf_ram / 2, std::max((std::uint64_t)4,
          (io_buf_ram / 2) / ((std::uint64_t)2 << 20)));

    // Initialize the reader of SA.
    typedef async_stream_reader<text_offset_type> sa_reader_type;
    sa_reader_type *sa_reader = new sa_reader_type(sa_filename,
        io_buf_ram / 2, std::max((std::uint64_t)4,
          (io_buf_ram / 2) / ((std::uint64_t)2 << 20)));

    // Compute final BWT.
    for (std::uint64_t j = 0; j < text_length; ++j) {
      std::uint64_t sa_j = sa_reader->read();
      if (sa_j == 0) bwt_writer->write((char_type)0);
      else {
        --sa_j;
        std::uint64_t block_id = sa_j / max_block_size;
        char_type bwt_j =
          bwt_subseq_reader->read_from_ith_file(block_id);
        bwt_writer->write(bwt_j);
      }
    }

    // Stop I/O threads.
    sa_reader->stop_reading();
    bwt_subseq_reader->stop_reading();

    // Update I/O volume.
    io_volume +=
      sa_reader->bytes_read() +
      bwt_subseq_reader->bytes_read() +
      bwt_writer->bytes_written();
    total_io_volume += io_volume;

    // Clean up.
    delete sa_reader;
    delete bwt_writer;
    delete bwt_subseq_reader;
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id)
      utils::file_delete(bwt_subsequences_filenames[block_id]);

    // Print summary.
    long double merge_bwt_subseq_time =
      utils::wclock() - merge_bwt_subseq_start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, "
        "total I/O = %.2Lfbytes/symbol\n",
        merge_bwt_subseq_time,
        ((1.L * io_volume) / (1L << 20)) / merge_bwt_subseq_time,
        (1.L * total_io_volume) / text_length);
  }

  // Clean up.
  delete[] sa_subsequences_filenames;
  delete[] bwt_subsequences_filenames;

  // Print summary.
  long double total_time = utils::wclock() - global_start;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  Absolute time = %.2Lfs\n", total_time);
  fprintf(stderr, "  Relative time = %.2Lfus/symbol\n",
      (1000000.0 * total_time) / text_length);
  fprintf(stderr, "  I/O volume = %lu bytes (%.2Lfbytes/symbol)\n",
      total_io_volume, (1.L * total_io_volume) / text_length);

#ifdef MONITOR_DISK_USAGE
  fprintf(stderr, "  Internal I/O volume counter = %lu\n",
      utils::get_current_io_volume());
#endif

  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));

#ifdef MONITOR_DISK_USAGE
  fprintf(stderr, "  Disk allocation: cur = %.2LfGiB, peak = %.2LfGiB\n",
      (1.L * utils::get_current_disk_allocation()) / (1UL << 30),
      (1.L * utils::get_peak_disk_allocation()) / (1UL << 30));
#endif

}

}  // namespace pem_bwt_private

template<typename char_type,
  typename text_offset_type>
void compute_bwt(
    std::string text_filename,
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t ram_use) {
  pem_bwt_private::compute_bwt<char_type, text_offset_type>(
      text_filename, sa_filename, output_filename, ram_use);
}

#endif  // __SRC_PEM_BWT_SRC_COMPUTE_BWT_HPP_INCLUDED
