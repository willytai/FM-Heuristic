#ifndef DLIST_H
#define DLIST_H

#include <cassert>
#include <iostream>
#include <climits>
#include <vector>
#include <algorithm>

using namespace std;

enum Group {  A = 0, B = 1, X = 2};

struct Cell {
    Cell(int id, Group g) : _next(NULL), _prev(NULL), _gain(0), _group(X), _locked(false) { _ID = id; _group = g; }
    ~Cell() {}

    void change_group() { if (_group == A) _group = B; else _group = A; }
    void set_gain(int g) { _gain = g; }
    void lock() { _locked = true; }
    void unlock() { _locked = false; }
    void increment_gain() { ++_gain; }
    void decrement_gain() { --_gain; }

    Cell*  _next;
    Cell*  _prev;
    int    _ID;
    int    _gain;
    Group  _group;
    bool   _locked;
};

class Dlist {
public:
    Dlist() : _size(0) { _head = NULL; }
    ~Dlist() {}

    void insert(Cell* c) {
        if (!_head) { this->init(c); return; }
        c->_prev = NULL;
        c->_next = _head;
        _head->_prev = c;
        _head = c;
        ++_size;
        _pick_candidate = _head;
    }

    void remove(Cell* c) {
#ifdef DEBUG_MODE
        if (!this->check_exists(c)) {
            cerr << "Cell " << c->_ID << " does not exist in list!!!" << endl;
            assert(0);
        }
#endif
        if (_size == 1) { this->clear(); return; }
        if (c == _head) _head = c->_next;
        if (c->_next) c->_next->_prev = c->_prev;
        if (c->_prev) c->_prev->_next = c->_next;
        --_size;
        _pick_candidate = _head;
    }

    void init(Cell* c) {
        _head = c;
        _head->_prev = _head->_next = NULL;
        _size = 1;
        _pick_candidate = _head;
    }

    void clear() {
        _head = NULL;
        _pick_candidate = NULL;
        _size = 0;
    }

    void print() const {
        Cell* tmp = _head;
        while (tmp) {
            cout << tmp->_ID << ' ';
            tmp = tmp->_next;
        } cout << endl;
    }

    void reset_pick_candidate() {
        _pick_candidate = _head;
    }

#ifdef DEBUG_MODE
    bool check_exists(Cell* c) {
        Cell* tmp = _head;
        while (tmp) {
            if (tmp == c) return true;
            tmp = tmp->_next;
        }
        return false;
    }
#endif

    int Gain() const { if (!_head) return INT_MIN; return _head->_gain; }
    int size() const { return _size; }

    Cell* pick() {
        Cell* tmp = _pick_candidate;
        _pick_candidate = _pick_candidate ? _pick_candidate->_next: _pick_candidate;
        return tmp;
        /* Cell* tmp = _head;
        while (id && tmp) {
            --id;
            tmp = tmp->_next;
        }
        return tmp;*/
    }

    Group group() const { return _head->_group; }

private:
    Cell*   _head;
    Cell*   _pick_candidate;
    int     _size;
};

class Bucket{
public:
    Bucket() : _size(0) {}
    ~Bucket() {}

    void set_gain_limit_and_resize(int limit) { _Pmax = limit; _data.resize(2*limit+1); }
    void insert(int gain, Cell* c) { (*this)[gain].insert(c); ++_size; }
    void remove(int gain, Cell* c) { (*this)[gain].remove(c); --_size; }
    void insert_locked(Cell* c) { ++_size; _locked_cells.push_back(c); }

    void clear_locked_cells() {
        for (auto it = _locked_cells.begin(); it != _locked_cells.end(); ++it) (*it)->unlock();
        _size -= _locked_cells.size();
        _locked_cells.clear();
    }

    void dump(ostream& os) {
        ::sort(_locked_cells.begin(), _locked_cells.end());
        for (auto it = _locked_cells.begin(); it != _locked_cells.end(); ++it) {
            os << 'c' << (*it)->_ID << ' ';
        }
        os << ';' << endl;
    }


    size_t size() const { return _size; }

    Dlist& operator[] (int gain) { assert(gain >= -_Pmax && gain <= _Pmax); return _data[gain+_Pmax]; }
    Dlist operator[] (int gain) const { assert(gain >= -_Pmax && gain <= _Pmax); return _data[gain+_Pmax]; }

    // debug functions
    void print() const {
        cerr << "[Bucket Content]" << endl;
        for (int i = _Pmax; i >= -_Pmax; --i) {
            cerr << "[gain " << i << "] ";
            (*this)[i].print();
        }
        cerr << endl;
    }

private:
    vector<Dlist>   _data;
    vector<Cell*>   _locked_cells;
    int             _Pmax;
    size_t          _size;
};

#endif /* DLIST_H */
