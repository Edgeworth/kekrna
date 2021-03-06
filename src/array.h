// Copyright 2016, Eliot Courtney.
//
// This file is part of kekrna.
//
// kekrna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// kekrna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with kekrna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef KEKRNA_ARRAY_H
#define KEKRNA_ARRAY_H

#include "common.h"

namespace kekrna {

template <typename T, unsigned int K>
struct array3d_t {
  typedef T ArrayType[K];

public:
  array3d_t() : data(nullptr), size(0) {}
  array3d_t(std::size_t size_, uint8_t init_val = MAX_E & 0xFF)
      : data(new T[size_ * size_ * K]), size(size_) {
    memset(data, init_val, sizeof(data[0]) * size_ * size_ * K);
  }
  ~array3d_t() { delete[] data; }

  array3d_t(const array3d_t&) = delete;
  array3d_t& operator=(const array3d_t&) = delete;
  array3d_t(array3d_t&& o) : data(nullptr), size(0) { *this = std::move(o); }

  array3d_t& operator=(array3d_t&& o) {
    delete[] data;
    data = o.data;
    size = o.size;
    o.data = nullptr;
    o.size = 0;
    return *this;
  }

  ArrayType* operator[](std::size_t idx) {
    return reinterpret_cast<ArrayType*>(&data[idx * size * K]);
  }

  const ArrayType* operator[](std::size_t idx) const {
    return reinterpret_cast<const ArrayType*>(&data[idx * size * K]);
  }

private:
  T* data;
  std::size_t size;
};

template <typename T, unsigned int K>
struct array2d_t {
public:
  array2d_t() : data(nullptr), size(0) {}
  ~array2d_t() { delete[] data; }

  array2d_t(std::size_t size_, uint8_t init_val = MAX_E & 0xFF)
      : data(new T[size_ * K]), size(size_) {
    memset(data, init_val, sizeof(data[0]) * size_ * K);
  }

  array2d_t(const array2d_t&) = delete;
  array2d_t& operator=(const array2d_t&) = delete;
  array2d_t(array2d_t&& o) : data(nullptr), size(0) { *this = std::move(o); }

  array2d_t& operator=(array2d_t&& o) {
    delete[] data;
    data = o.data;
    size = o.size;
    o.data = nullptr;
    o.size = 0;
    return *this;
  }

  T* operator[](std::size_t idx) { return &data[idx * K]; }

  const T* operator[](std::size_t idx) const { return &data[idx * K]; }

private:
  T* data;
  std::size_t size;
};
}

#endif  // KEKRNA_ARRAY_H
