#include "Solver.h"
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

void Solver::read(char* filename) {
    cerr << "[constructing Cell array and Net array from " << filename << "]" << endl;
    ifstream file(filename);
    string buff;
    bool A_check = false;
    bool B_check = false;
    while (getline(file, buff)) {
        if (!buff.length()) continue;
        if (!_balance_degree) {
            _balance_degree = stof(buff);
            cerr << "\t> setting balance degree to " << _balance_degree << endl;
            continue;
        }
        istringstream ss(buff);
        string temp;
        int curNetID;
        int curCellID;
        while (ss >> temp) {
            if (temp == ";") {
                if (A_check && B_check) ++_cutsize;
                break;
            }
            if (temp == "NET") continue;
            if (temp[0] == 'n') {
                curNetID = stoi(temp.substr(1,temp.size()-1));
                A_check = B_check = false;
                continue;
            }
            if (temp[0] == 'c') curCellID = stoi(temp.substr(1,temp.size()-1));

            if (_cell_array.size() < size_t(curCellID + 1)) _cell_array.resize(curCellID+1);
            if (_net_array.size() < size_t(curNetID + 1)) _net_array.resize(curNetID+1);

            if (this->net_in_cell_array(curCellID, curNetID)) continue;

            _cell_array[curCellID].push_back(curNetID);
            _net_array[curNetID].push_back(curCellID);

            if (_numCell < curCellID) _numCell = curCellID;
            if (_numNet < curNetID) _numNet = curNetID;

            // create cell pointers here
            // maybe some ID's are not used
            // initialize partition now
            if (_cell_ptr.size() < _numCell+1) _cell_ptr.resize(_numCell+1, NULL);
            if (curCellID&1) {
                _cell_ptr[curCellID] = new Cell(curCellID, A);
                A_check = true;
            }
            else {
                _cell_ptr[curCellID] = new Cell(curCellID, B);
                B_check = true;
            }
        }
    }

    // compute _Pmax
    for (auto it = _cell_array.begin(); it != _cell_array.end(); ++it) {
        if (_Pmax < (*it).size()) _Pmax = (*it).size();
    }
    
    cerr << "\t> Number of Cells: " << _numCell << endl;
    cerr << "\t> Number of Nets:  " << _numNet << endl;
    cerr << "\t> initial cutsize: " << _cutsize << endl;
}

void Solver::solve() {
    _cell_gain_pairs.clear();
    this->construct_balance_criterion();
    this->initBucketList();
    int iteration = 0;
    cerr << "[Solving]" << endl;
    while (true) {
        cout << "\t> iteration " << setw(3) << right << ++iteration << " -> gain ";
        while (true) {
            this->moveMaxGainCell();
            if (!this->update_gain()) break;
        }
        if (!this->compute_max_gain()) break;
    }
}

void Solver::construct_balance_criterion() {
    cerr << "[constructing balance criterion]" << endl;
    _min_limit = (1 - _balance_degree) / 2 * _numCell;
    _max_limit = (1 + _balance_degree) / 2 * _numCell;
    cerr << "\t> balance range [" << _min_limit << ", " << _max_limit << "]" << endl;
}

void Solver::initBucketList() {

    // resize bucketlist
    _Bucket.clear();
    _Bucket.set_gain_limit_and_resize(2*_Pmax+1);

    // resize net distribution and initialize
    _NetADistribution.resize(_numNet+1);
    _NetBDistribution.resize(_numNet+1);
    for (int netID = 1; netID <= _numNet; ++netID) {
        int A_count = 0;
        int B_count = 0;
        for (auto cell_it = _net_array[netID].begin(); cell_it != _net_array[netID].end(); ++cell_it) {
            int cellID = *cell_it;
            if (_cell_ptr[cellID]->_group == A) ++A_count;
            else ++B_count;
        }
        _NetADistribution[netID] = A_count;
        _NetBDistribution[netID] = B_count;
    }

    // init cell gain, place in bucket and set _maxGainPtr
    for (int cellID = 1; cellID <= _numCell; ++cellID) {
        int gain = 0;

        // clearify From Block
        Group F = _cell_ptr[cellID]->_group;

        for (auto net_it = _cell_array[cellID].begin(); net_it != _cell_array[cellID].end(); ++net_it) {

            int netID = *net_it;
            if (F == A) {
                if (_NetADistribution[netID] == 1) ++gain;
                if (_NetBDistribution[netID] == 0) --gain;
            }
            else {
                if (_NetBDistribution[netID] == 1) ++gain;
                if (_NetADistribution[netID] == 0) --gain;
            }
        }
        _cell_ptr[cellID]->set_gain(gain);

        // update bucketlist
        _Bucket.insert(gain, _cell_ptr[cellID]);
    }
    this->update_max_gain_pointer();
}

void Solver::moveMaxGainCell() {
    // pick a legal cell
    Cell* candidate = _maxGainPtr->pick();
    this->assign_basecell(candidate, candidate->_group);
    while (true) {
        if (!this->balance_checking(_BaseCell)) candidate = _maxGainPtr->pick();
        else break;
        if (!candidate) {
            this->move_max_gain_pointer();
            candidate = _maxGainPtr->pick();
        }
        this->assign_basecell(candidate, candidate->_group);
    }
    // remove it from bucket
    _Bucket.remove(_BaseCell->_gain, _BaseCell);
    // lock and change group
    _BaseCell->lock();
    _BaseCell->change_group();
    // put into locked group
    _Bucket.insert_locked(_BaseCell);
    // record
    _cell_gain_pairs.push_back(pair<Cell*,int>(_BaseCell, _BaseCell->_gain));
}

bool Solver::balance_checking(Cell* c) {
    int A_size = _Bucket.A_size();
    int B_size = _Bucket.B_size();
    if (c->_group == A) {
        --A_size;
        ++B_size;
    }
    else {
        ++A_size;
        --B_size;
    }
    if (A_size < _min_limit || A_size > _max_limit) return false;
    else return true;
}

bool Solver::update_gain() {
    vector<int>* F_distribution = &_NetADistribution;
    vector<int>* T_distribution = &_NetBDistribution;
    if (_BaseCell_F == B) ::swap(F_distribution, T_distribution);

    for (auto net_it = _cell_array[_BaseCell->_ID].begin(); net_it != _cell_array[_BaseCell->_ID].end(); ++net_it) {
        int netID = *net_it;
        for (auto cell_it = _net_array[netID].begin(); cell_it != _net_array[netID].end(); ++cell_it) {
            if (_cell_ptr[*cell_it]->_locked) continue;
            int cellID = *cell_it;
            // increment gains of all free cells on netID
            if (T_distribution->at(netID) == 0) {
                _Bucket.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                _cell_ptr[cellID]->increment_gain();
                _Bucket.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
            }
            // decrement gains of free cells in T on netID
            if (T_distribution->at(netID) == 1 && _cell_ptr[cellID]->_group != _BaseCell_F) {
                _Bucket.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                _cell_ptr[cellID]->decrement_gain();
                _Bucket.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
            }
        }
    }
    this->update_net_distribution(F_distribution, T_distribution);
    for (auto net_it = _cell_array[_BaseCell->_ID].begin(); net_it != _cell_array[_BaseCell->_ID].end(); ++net_it) {
        int netID = *net_it;
        for (auto cell_it = _net_array[netID].begin(); cell_it != _net_array[netID].end(); ++cell_it) {
            if (_cell_ptr[*cell_it]->_locked) continue;
            int cellID = *cell_it;
            // decrement gains of free cells in F on netID
            if (F_distribution->at(netID) == 0) {
                _Bucket.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                _cell_ptr[cellID]->decrement_gain();
                _Bucket.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
            }
            // increment gains of all free cells on netID
            if (F_distribution->at(netID) == 1 && _cell_ptr[cellID]->_group == _BaseCell_F) {
                _Bucket.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                _cell_ptr[cellID]->increment_gain();
                _Bucket.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
            }
        }
    }
    return this->update_max_gain_pointer();
}

bool Solver::update_max_gain_pointer() {
    for (int gain = _Pmax; gain >= -_Pmax; --gain) {
        if (_Bucket[gain].size()) {
            _maxGainPtr = &_Bucket[gain];
            _maxGainPtr->reset_pick_candidate();
            return true;
        }
    }
    return false;
}

bool Solver::compute_max_gain() {
    int gain = 0;
    int max_gain = 0;
    int k = 0;
    for (unsigned int i = 0; i < _cell_gain_pairs.size(); ++i) {
        gain += _cell_gain_pairs[i].second;
        if (max_gain < gain) {
            k = i;
            max_gain = gain;
        }
    }
    if (max_gain > 0) {
        this->update_cutsize(max_gain);
        this->apply_change(k);
        cout << max_gain << endl;
        return true;
    }
    else {
        cout << "0" << endl;
        return false;
    }
}

void Solver::update_cutsize(const int& gain) {
    _cutsize -= gain;
}

void Solver::apply_change(int k) {
    // unlock all the cells
    for (int i = 1; i <= k; ++i) {
        _cell_gain_pairs[i].first->unlock();
    }
    for (unsigned int i = k+1; i < _cell_gain_pairs.size(); ++i) {
        _cell_gain_pairs[i].first->unlock();
        _cell_gain_pairs[i].first->change_group();
    }
    _cell_gain_pairs.clear();
    this->initBucketList();
}

void Solver::dump(ostream& os) {
    cerr << "[Dumping result]" << endl;
    os << "Cutsize = " << _cutsize << endl;
    os << "G1 " << _Bucket.A_size() << endl;
    _Bucket.dumpA(os);
    os << "G2 " << _Bucket.B_size() << endl;
    _Bucket.dumpB(os);
}

inline void Solver::move_max_gain_pointer() {
    int curGain = _maxGainPtr->Gain();
    while (!_Bucket[--curGain].size()) {}
    _maxGainPtr = &_Bucket[curGain];
    _maxGainPtr->reset_pick_candidate();
}

inline void Solver::update_net_distribution(vector<int>* F, vector<int>* T) {
    // update distribution for every net on BaseCell
    for (auto net_it = _cell_array[_BaseCell->_ID].begin(); net_it != _cell_array[_BaseCell->_ID].end(); ++net_it) {
        --(F->at(*net_it));
        ++(T->at(*net_it));
    }
}

inline void Solver::assign_basecell(Cell* c, Group g) {
    _BaseCell = c;
    _BaseCell_F = g;
}

bool Solver::net_in_cell_array(int curCellID, int curNetID) {
    for (auto it = _net_array[curNetID].begin(); it != _net_array[curNetID].end(); ++it) {
        if (*it == curCellID) return true;
    }
    return false;
}

#ifdef DEBUG_MODE
void Solver::debug_net_dist() const {
    const vector<int>* F = &_NetADistribution;
    const vector<int>* T = &_NetBDistribution;
    if (_BaseCell_F == B) ::swap(F, T);
    cout << "    F | T |" << endl;
    for (int netID = 1; netID <= _numNet; ++netID) {
        cout << 'n' << netID << "  " << F->at(netID) << " | " << T->at(netID) << " | " << endl;
    }
}
#endif
