// sqlite3ppext.cpp
//
// The MIT License
//
// Copyright (c) 2015 Wongoo Lee (iwongu at gmail dot com)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include <cstring>

#include "sqlite3ppext.h"

namespace sqlite3pp
{
  namespace ext
  {

    namespace
    {


    } // namespace


    context::context(sqlite3_context* ctx, int nargs, sqlite3_value** values)
      : ctx_(ctx), nargs_(nargs), values_(values)
    {
    }


    int context::get(int idx, int) const
    {
      return sqlite3_value_int(values_[idx]);
    }

    double context::get(int idx, double) const
    {
      return sqlite3_value_double(values_[idx]);
    }

    long long int context::get(int idx, long long int) const
    {
      return sqlite3_value_int64(values_[idx]);
    }

    char const* context::get(int idx, char const*) const
    {
      return reinterpret_cast<char const*>(sqlite3_value_text(values_[idx]));
    }

    std::string context::get(int idx, std::string) const
    {
      return get(idx, (char const*)0);
    }

    void const* context::get(int idx, void const*) const
    {
      return sqlite3_value_blob(values_[idx]);
    }



    void context::result(int value)
    {
      sqlite3_result_int(ctx_, value);
    }

    void context::result(double value)
    {
      sqlite3_result_double(ctx_, value);
    }

    void context::result(long long int value)
    {
      sqlite3_result_int64(ctx_, value);
    }

    void context::result(std::string const& value)
    {
      result(value.c_str(), false);
    }

    void context::result(char const* value, bool fstatic)
    {
      sqlite3_result_text(ctx_, value, std::strlen(value), fstatic ? SQLITE_STATIC : SQLITE_TRANSIENT);
    }

    void context::result(void const* value, int n, bool fstatic)
    {
      sqlite3_result_blob(ctx_, value, n, fstatic ? SQLITE_STATIC : SQLITE_TRANSIENT);
    }

    void context::result()
    {
      sqlite3_result_null(ctx_);
    }

    void context::result(null_type)
    {
      sqlite3_result_null(ctx_);
    }

    void* context::aggregate_data(int size)
    {
      return sqlite3_aggregate_context(ctx_, size);
    }

    int context::aggregate_count()
    {
      return sqlite3_aggregate_count(ctx_);
    }

    function::function(database& db) : db_(db.db_)
    {
    }



    aggregate::aggregate(database& db) : db_(db.db_)
    {
    }


  } // namespace ext

} // namespace sqlite3pp
