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
    void set_group(Group g) { _group = g; }
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
    Dlist() : _size(0), _gain(INT_MAX) { _head = NULL; }
    ~Dlist() {}

    void setGain(int g) { _gain = g; }

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

    int Gain() const { return _gain; }
    int size() const { return _size; }

    Cell* pick() {
        Cell* tmp = _pick_candidate;
        _pick_candidate = _pick_candidate ? _pick_candidate->_next: _pick_candidate;
        return tmp;
    }

    Group group() const { return _head->_group; }

private:
    Cell*   _head;
    Cell*   _pick_candidate;
    int     _size;
    int     _gain;
};

class Bucket{
public:
    Bucket() : _size(0), _A_size(0), _B_size(0) {}
    ~Bucket() {}

    void insert(int gain, Cell* c) { (*this)[gain].insert(c); ++_size; if (c->_group == A) ++_A_size; else ++_B_size; }
    void remove(int gain, Cell* c) { (*this)[gain].remove(c); --_size; if (c->_group == A) --_A_size; else --_B_size; }
    void insert_locked(Cell* c) { ++_size; _locked_cells.push_back(c); if (c->_group == A) ++_A_size; else ++_B_size; }

    void set_gain_limit_and_resize(int limit) {
        _Pmax = limit;
        _data.resize(2*limit+1);
        for (unsigned int i = 0; i < _data.size(); ++i) _data[i].setGain(i-limit);
    }

    void clear_locked_cells() {
        for (auto it = _locked_cells.begin(); it != _locked_cells.end(); ++it) (*it)->unlock();
        _size -= _locked_cells.size();
        _locked_cells.clear();
    }

    void dumpA(ostream& os) {
        ::sort(_locked_cells.begin(), _locked_cells.end());
        for (auto it = _locked_cells.begin(); it != _locked_cells.end(); ++it) {
            if ((*it)->_group == B) continue;
            os << 'c' << (*it)->_ID << ' ';
        }
        os << ';' << endl;
    }

    void dumpB(ostream& os) {
        for (auto it = _locked_cells.begin(); it != _locked_cells.end(); ++it) {
            if ((*it)->_group == A) continue;
            os << 'c' << (*it)->_ID << ' ';
        }
        os << ';' << endl;
    }

    void clear() {
        _data.clear();
        _locked_cells.clear();
        _size = _A_size = _B_size = 0;
    }

    bool no_more_free_cells() const {
        return _locked_cells.size() == _size;
    }

    size_t size()   const { return _size; }
    size_t A_size() const { return _A_size; }
    size_t B_size() const { return _B_size; }

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

    void lock_insert(Cell* c) {
        _locked_cells.push_back(c);
        if (c->_group == A) _A_size++;
        else _B_size++;
    }

private:
    vector<Dlist>   _data;
    vector<Cell*>   _locked_cells;
    int             _Pmax;
    size_t          _size;
    size_t          _A_size;
    size_t          _B_size;
};

#endif /* DLIST_H */
