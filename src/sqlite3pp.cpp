// sqlite3pp.cpp
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
#include <memory>

#include "sqlite3pp.h"

namespace sqlite3pp{

  null_type ignore;
  namespace{} // namespace
  database::database(char const* dbname, int flags, char const* vfs) : db_(nullptr) {
    if (dbname) {
      auto rc = connect(dbname, flags, vfs);
      if (rc != SQLITE_OK)
        throw database_error("can't connect database");
    }
  }

  database::~database()
  {
    disconnect();
  }

  int database::connect(char const* dbname, int flags, char const* vfs)
  {
    disconnect();

    return sqlite3_open_v2(dbname, &db_, flags, vfs);
  }

  int database::disconnect()
  {
    auto rc = SQLITE_OK;
    if (db_) {
      rc = sqlite3_close(db_);
      db_ = nullptr;
    }

    return rc;
  }




  int database::execute(char const* sql)
  {
    return sqlite3_exec(db_, sql, 0, 0, 0);
  }


  statement::statement(database& db, char const* stmt) : db_(db), stmt_(0), tail_(0)
  {
    if (stmt) {
      auto rc = prepare(stmt);
      if (rc != SQLITE_OK)
        throw database_error(db_);
    }
  }

  statement::~statement()
  {
    auto rc = finish();
    if (rc != SQLITE_OK)
      throw database_error(db_);
  }

  int statement::prepare(char const* stmt)
  {
    auto rc = finish();
    if (rc != SQLITE_OK)
      return rc;

    return prepare_impl(stmt);
  }

  int statement::prepare_impl(char const* stmt)
  {
    return sqlite3_prepare(db_.db_, stmt, std::strlen(stmt), &stmt_, &tail_);
  }

  int statement::finish()
  {
    auto rc = SQLITE_OK;
    if (stmt_) {
      rc = finish_impl(stmt_);
      stmt_ = nullptr;
    }
    tail_ = nullptr;

    return rc;
  }

  int statement::finish_impl(sqlite3_stmt* stmt)
  {
    return sqlite3_finalize(stmt);
  }

  int statement::step()
  {
    return sqlite3_step(stmt_);
  }


  int statement::bind(int idx, int value)
  {
    return sqlite3_bind_int(stmt_, idx, value);
  }

  int statement::bind(int idx, double value)
  {
    return sqlite3_bind_double(stmt_, idx, value);
  }

  int statement::bind(int idx, long long int value)
  {
    return sqlite3_bind_int64(stmt_, idx, value);
  }

  int statement::bind(int idx, char const* value, bool fstatic)
  {
    return sqlite3_bind_text(stmt_, idx, value, std::strlen(value), fstatic ? SQLITE_STATIC : SQLITE_TRANSIENT);
  }

  int statement::bind(int idx, void const* value, int n, bool fstatic)
  {
    return sqlite3_bind_blob(stmt_, idx, value, n, fstatic ? SQLITE_STATIC : SQLITE_TRANSIENT);
  }

  int statement::bind(int idx)
  {
    return sqlite3_bind_null(stmt_, idx);
  }

  int statement::bind(int idx, null_type)
  {
    return bind(idx);
  }

  int statement::bind(char const* name, int value)
  {
    auto idx = sqlite3_bind_parameter_index(stmt_, name);
    return bind(idx, value);
  }

  int statement::bind(char const* name, double value)
  {
    auto idx = sqlite3_bind_parameter_index(stmt_, name);
    return bind(idx, value);
  }

  int statement::bind(char const* name, long long int value)
  {
    auto idx = sqlite3_bind_parameter_index(stmt_, name);
    return bind(idx, value);
  }

  int statement::bind(char const* name, char const* value, bool fstatic)
  {
    auto idx = sqlite3_bind_parameter_index(stmt_, name);
    return bind(idx, value, fstatic);
  }

  int statement::bind(char const* name, void const* value, int n, bool fstatic)
  {
    auto idx = sqlite3_bind_parameter_index(stmt_, name);
    return bind(idx, value, n, fstatic);
  }

  int statement::bind(char const* name)
  {
    auto idx = sqlite3_bind_parameter_index(stmt_, name);
    return bind(idx);
  }

  int statement::bind(char const* name, null_type)
  {
    return bind(name);
  }


  command::bindstream::bindstream(command& cmd, int idx) : cmd_(cmd), idx_(idx)
  {
  }

  command::command(database& db, char const* stmt) : statement(db, stmt)
  {
  }

  command::bindstream command::binder(int idx)
  {
    return bindstream(*this, idx);
  }

  int command::execute()
  {
    auto rc = step();
    if (rc == SQLITE_DONE) rc = SQLITE_OK;

    return rc;
  }


  query::rows::getstream::getstream(rows* rws, int idx) : rws_(rws), idx_(idx)
  {
  }

  query::rows::rows(sqlite3_stmt* stmt) : stmt_(stmt)
  {
  }

  int query::rows::data_count() const
  {
    return sqlite3_data_count(stmt_);
  }

  int query::rows::column_type(int idx) const
  {
    return sqlite3_column_type(stmt_, idx);
  }

  int query::rows::column_bytes(int idx) const
  {
    return sqlite3_column_bytes(stmt_, idx);
  }

  int query::rows::get(int idx, int) const
  {
    return sqlite3_column_int(stmt_, idx);
  }

  double query::rows::get(int idx, double) const
  {
    return sqlite3_column_double(stmt_, idx);
  }

  long long int query::rows::get(int idx, long long int) const
  {
    return sqlite3_column_int64(stmt_, idx);
  }

  char const* query::rows::get(int idx, char const*) const
  {
    return reinterpret_cast<char const*>(sqlite3_column_text(stmt_, idx));
  }

  std::string query::rows::get(int idx, std::string) const
  {
    return get(idx, (char const*)0);
  }

  void const* query::rows::get(int idx, void const*) const
  {
    return sqlite3_column_blob(stmt_, idx);
  }

  null_type query::rows::get(int idx, null_type) const
  {
    return ignore;
  }
  query::rows::getstream query::rows::getter(int idx)
  {
    return getstream(this, idx);
  }

  query::query_iterator::query_iterator() : cmd_(0)
  {
    rc_ = SQLITE_DONE;
  }

  query::query_iterator::query_iterator(query* cmd) : cmd_(cmd)
  {
    rc_ = cmd_->step();
    if (rc_ != SQLITE_ROW && rc_ != SQLITE_DONE)
      throw database_error(cmd_->db_);
  }

  bool query::query_iterator::operator==(query::query_iterator const& other) const
  {
    return rc_ == other.rc_;
  }

  bool query::query_iterator::operator!=(query::query_iterator const& other) const
  {
    return rc_ != other.rc_;
  }

  query::query_iterator& query::query_iterator::operator++()
  {
    rc_ = cmd_->step();
    if (rc_ != SQLITE_ROW && rc_ != SQLITE_DONE)
      throw database_error(cmd_->db_);
    return *this;
  }

  query::query_iterator::value_type query::query_iterator::operator*() const
  {
    return rows(cmd_->stmt_);
  }

  query::query(database& db, char const* stmt) : statement(db, stmt)
  {
  }





  query::iterator query::begin()
  {
    return query_iterator(this);
  }

  query::iterator query::end()
  {
    return query_iterator();
  }
  transaction::transaction(database& db, bool fcommit, bool freserve) : db_(&db), fcommit_(fcommit){
    db_->execute(freserve ? "BEGIN IMMEDIATE" : "BEGIN");
  }
  transaction::~transaction()
  {
    if (db_) {
      auto rc = db_->execute(fcommit_ ? "COMMIT" : "ROLLBACK");
      if (rc != SQLITE_OK)
	throw database_error(*db_);
    }
  }
  database_error::database_error(char const* msg) : std::runtime_error(msg) { }
  database_error::database_error(database& db) : std::runtime_error(sqlite3_errmsg(db.db_)) { }
} // namespace sqlite3pp
