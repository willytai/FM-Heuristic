#include "Solver.h"
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

void Solver::read(char* filename) {
    cerr << "[constructing Cell array and Net array]" << endl;
    ifstream file(filename);
    string buff;
    while (getline(file, buff)) {
        if (!buff.length()) break;
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
            if (temp == ";") break;
            if (temp == "NET") continue;
            if (temp[0] == 'n') {
                curNetID = stoi(temp.substr(1,temp.size()-1));
                continue;
            }
            if (temp[0] == 'c') curCellID = stoi(temp.substr(1,temp.size()-1));
            if (_cell_array.size() < size_t(curCellID + 1)) _cell_array.resize(curCellID+1);
            if (_net_array.size() < size_t(curNetID + 1)) _net_array.resize(curNetID+1);
            _cell_array[curCellID].push_back(curNetID);
            _net_array[curNetID].push_back(curCellID);

            if (_numCell < curCellID) _numCell = curCellID;
            if (_numNet < curNetID) _numNet = curNetID;
        }
    }

    // compute _Pmax
    cerr << "\t> computing Pmax";
    for (auto it = _cell_array.begin(); it != _cell_array.end(); ++it) {
        if (_Pmax < (*it).size()) _Pmax = (*it).size();
    }
    cerr << " " << _Pmax << endl;
}

void Solver::solve() {
    _gain_history.clear();
    this->construct_balance_criterion();
    this->initPartition();
    this->initBucketList();
    int iteration = 0;
    while (true) {
        cerr << "[main process, iteration " << ++iteration << "]" << endl;
        _cell_gain_pairs.clear();
        int step = 0;
        while (true) {
            cerr << "\r\t> step " << ++step << flush;
            this->moveMaxGainCell();
            if (!this->update_gain()) break;
        } cerr << endl;
        if (!this->compute_max_gain()) break;
    }
    this->print_hisotry();
}

void Solver::construct_balance_criterion() {
    cerr << "[constructing balance criterion]" << endl;
    _min_limit = (1 - _balance_degree) / 2 * _numCell;
    _max_limit = (1 + _balance_degree) / 2 * _numCell;
    cerr << "\t> balance range [" << _min_limit << ", " << _max_limit << "]" << endl;
}

void Solver::initPartition() {
    cerr << "[initialzing partition and computing original cutsize]" << endl;
    int half = _numCell / 2;
    _cell_ptr.resize(_numCell+1, NULL);
    cerr << "\t> creating cell pointers for group A" << endl;
    for (int i = 1; i <= half; ++i) {
        _cell_ptr[i] = new Cell(i, A);
    }
    cerr << "\t> creating cell pointers for group B" << endl;
    for (int i = half+1; i <= _numCell; ++i) {
        _cell_ptr[i] = new Cell(i, B);
    }
    vector<bool> checker;
    checker.resize(_numNet+1, false);
    for (int i = 1; i <= half; ++i) {
        for (auto it = _cell_array[i].begin(); it != _cell_array[i].end(); ++it) {
            checker[*it] = true;
        }
    }
    for (int i = half+1; i <= _numCell; ++i) {
        for (auto it = _cell_array[i].begin(); it != _cell_array[i].end(); ++it) {
            if (checker[*it]) {
                _cutsize += 1;
                checker[*it] = false;
            }
        }
    }
    cerr << "\t> initial cutsize: " << _cutsize << endl;
}

void Solver::initBucketList() {
    cerr << endl << "[initialzing bucketlist for group A and B]" << endl;

    // resize bucketlist
    _Bucket_A.set_gain_limit_and_resize(2*_Pmax+1);
    _Bucket_B.set_gain_limit_and_resize(2*_Pmax+1);
    _Bucket_A.clear_locked_cells();
    _Bucket_B.clear_locked_cells();

    // resize net distribution
    _NetADistribution.resize(_numNet+1);
    _NetBDistribution.resize(_numNet+1);

    // init cell gain, place in the corresponding bucket and set _maxGainPtr
    _maxGainPtr = &_Bucket_A[-_Pmax];
    for (int cellID = 1; cellID <= _numCell; ++cellID) {
        cerr << "\r\t> computing gain for cell " << cellID << '/' << _numCell << flush;
        int gain = 0;

        // clearify From Block
        Group F = _cell_ptr[cellID]->_group;

        for (auto net_it = _cell_array[cellID].begin(); net_it != _cell_array[cellID].end(); ++net_it) {

            // colculate number of cells in F and T respectively
            int netID   = *net_it;
            int F_count = 0;
            int T_count = 0;
            for (auto cell_it = _net_array[netID].begin(); cell_it != _net_array[netID].end(); ++cell_it) {
                if (_cell_ptr[*cell_it]->_group == F) ++F_count;
                else ++T_count;
            }
            if (F_count == 1) ++gain;
            if (T_count == 0) --gain;

            // init net distribution
            _NetADistribution[netID] = F_count;
            _NetBDistribution[netID] = T_count;
            if (F == B) ::swap(_NetADistribution[netID], _NetBDistribution[netID]);
        }
        _cell_ptr[cellID]->set_gain(gain);

        // update bucketlist and _maxGainPtr
        if (F == A) {
            _Bucket_A.insert(gain, _cell_ptr[cellID]);
            if (_maxGainPtr->Gain() < gain) {
                _maxGainPtr = &_Bucket_A[gain];
            }
        }
        else {
            _Bucket_B.insert(gain, _cell_ptr[cellID]);
            if (_maxGainPtr->Gain() < gain) {
                _maxGainPtr = &_Bucket_B[gain];
            }
        }
    } cerr << endl;
    cerr << "\t> Initial Distribution: A " << _Bucket_A.size() << ", B " << _Bucket_B.size() << endl;
}

void Solver::moveMaxGainCell() {
    // pick a legal cell
    cerr << "\tsearching for BaseCell...";
    this->assign_basecell(_maxGainPtr->pick(), _maxGainPtr->group()); assert(BaseCell);
    Group ref = BaseCell->_group;
    int times = 0;
    while (true) {
        ++times;
        if (!this->balance_checking(BaseCell)) {
            // cerr << "\tBaseCell : " << BaseCell->_ID << " with gain " << BaseCell->_gain << " cannot be swapped." << endl;
            this->assign_basecell(_maxGainPtr->pick(), _maxGainPtr->group());
        }
        else break;
        if (!BaseCell) {
            // cerr << "\tsearching for the second max gain" << endl;
            this->move_max_gain_pointer(ref);
            this->assign_basecell(_maxGainPtr->pick(), _maxGainPtr->group());
        }
    }
    // cerr << "done. " << times << endl;
    // cerr << "Swapping Cell " << BaseCell->_ID << " from group " << BaseCell->_group << " to group " << 1-int(BaseCell->_group) << endl;
    // remove it from bucket and put into locked group
    if (BaseCell->_group == A) {
        _Bucket_A.remove(BaseCell->_gain, BaseCell);
        _Bucket_B.insert_locked(BaseCell);
    }
    if (BaseCell->_group == B) {
        _Bucket_B.remove(BaseCell->_gain, BaseCell);
        _Bucket_A.insert_locked(BaseCell);
    }
    // lock and change group
    BaseCell->lock();
    BaseCell->change_group();
    // record
    _cell_gain_pairs.push_back(pair<Cell*,int>(BaseCell, BaseCell->_gain));
}

void Solver::move_max_gain_pointer(Group ref) {
    int curGain = _maxGainPtr->Gain();
    if (_maxGainPtr->group() == A && _Bucket_B[curGain].size() && _maxGainPtr->group() == ref) {
        _maxGainPtr = &_Bucket_B[curGain];
        _maxGainPtr->reset_pick_candidate();
        return;
    }
    if (_maxGainPtr->group() == B && _Bucket_A[curGain].size() && _maxGainPtr->group() == ref) {
        _maxGainPtr = &_Bucket_A[curGain];
        _maxGainPtr->reset_pick_candidate();
        return;
    }
    while (true) {
        --curGain;
        /* if (curGain < -_Pmax) {
            cerr << "something wrong!!!!!!" << endl;
            assert(0);
        }*/
        if (_Bucket_A[curGain].size()) {
            _maxGainPtr = &_Bucket_A[curGain];
            _maxGainPtr->reset_pick_candidate();
            return;
        }
        if (_Bucket_B[curGain].size()) {
            _maxGainPtr = &_Bucket_B[curGain];
            _maxGainPtr->reset_pick_candidate();
            return;
        }
    }
}

bool Solver::balance_checking(Cell* c) {
    // assert(c);
    // assert(!c->_locked);
    int A_size = _Bucket_A.size();
    int B_size = _Bucket_B.size();
    if (c->_group == A) {
        --A_size;
        ++B_size;
    }
    else {
        ++A_size;
        --B_size;
    }
    if (A_size < _min_limit || A_size > _max_limit) {
        // cerr << "balance checking failed. Cell " << c->_ID << " swapped to group " << 1-int(c->_group) << " makes the size of group 0 illegal (" << A_size << " is not in the range of [" << _min_limit << ", " << _max_limit << "]" << endl;
        return false;
    }
    else return true;
}

bool Solver::update_gain() {
    vector<int>* F_distribution = &_NetADistribution;
    vector<int>* T_distribution = &_NetBDistribution;
    if (_BaseCell_F == B) ::swap(F_distribution, T_distribution);

    for (auto net_it = _cell_array[BaseCell->_ID].begin(); net_it != _cell_array[BaseCell->_ID].end(); ++net_it) {
        int netID = *net_it;
        if (T_distribution->at(netID) == 0) { // increment gains of all free cells on netID
            for (auto cell_it = _net_array[netID].begin(); cell_it != _net_array[netID].end(); ++cell_it) {
                int cellID = *cell_it;
                if (_cell_ptr[cellID]->_locked) continue;
                else { // remember to update bucketlist
                    // cerr << "cell " << cellID << " incremented due to T(all)" << endl;
                    if (_cell_ptr[cellID]->_group == A) {
                        _Bucket_A.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                        _cell_ptr[cellID]->increment_gain();
                        _Bucket_A.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                    }
                    else {
                        _Bucket_B.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                        _cell_ptr[cellID]->increment_gain();
                        _Bucket_B.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                    }
                }
            }
        }
        else if (T_distribution->at(netID) == 1) { // decrement gains of the only T cell on netID
            for (auto cell_it = _net_array[netID].begin(); cell_it != _net_array[netID].end(); ++cell_it) {
                int cellID = *cell_it;
                if (_cell_ptr[cellID]->_locked || _cell_ptr[cellID]->_group == _BaseCell_F) continue;
                else { // remember to update bucketlist
                    // cerr << "cell " << cellID << " decremented due to T(only T)" << endl;
                    if (_cell_ptr[cellID]->_group == A) {
                        _Bucket_A.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                        _cell_ptr[cellID]->decrement_gain();
                        _Bucket_A.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                    }
                    else {
                        _Bucket_B.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                        _cell_ptr[cellID]->decrement_gain();
                        _Bucket_B.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                    }
                }
            }
        }
    }
    // cerr << "before" << endl;
    // this->debug_net_dist();
    this->update_net_distribution(F_distribution, T_distribution);
    // cerr << "after" << endl;
    // this->debug_net_dist();
    // assert(0);
    for (auto net_it = _cell_array[BaseCell->_ID].begin(); net_it != _cell_array[BaseCell->_ID].end(); ++net_it) {
        int netID = *net_it;
        if (F_distribution->at(netID) == 0) { // decrement gains of all free cells on netID
            for (auto cell_it = _net_array[netID].begin(); cell_it != _net_array[netID].end(); ++cell_it) {
                int cellID = *cell_it;
                if (_cell_ptr[cellID]->_locked) continue;
                else { // remember to update bucketlist
                    // cerr << "cell " << cellID << " decremented due to F'(all)" << endl;
                    if (_cell_ptr[cellID]->_group == A) {
                        _Bucket_A.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                        _cell_ptr[cellID]->decrement_gain();
                        _Bucket_A.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                    }
                    else {
                        _Bucket_B.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                        _cell_ptr[cellID]->decrement_gain();
                        _Bucket_B.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                    }
                }
            }
        }
        else if (F_distribution->at(netID) == 1) { // increment gains of the only F cell on netID
            for (auto cell_it = _net_array[netID].begin(); cell_it != _net_array[netID].end(); ++cell_it) {
                int cellID = *cell_it;
                if (_cell_ptr[cellID]->_locked || _cell_ptr[cellID]->_group != _BaseCell_F) continue;
                else { // remember to update bucketlist
                    // cerr << "cell " << cellID << " incremented due to F'(only F)" << endl;
                    if (_cell_ptr[cellID]->_group == A) {
                        _Bucket_A.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                        _cell_ptr[cellID]->increment_gain();
                        _Bucket_A.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                    }
                    else {
                        _Bucket_B.remove(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                        _cell_ptr[cellID]->increment_gain();
                        _Bucket_B.insert(_cell_ptr[cellID]->_gain, _cell_ptr[cellID]);
                    }
                }
            }
        }
    }
    // _Bucket_A.print();
    // _Bucket_B.print();
    // assert(0);
    return this->update_max_gain_pointer();
}

void Solver::update_net_distribution(vector<int>* F, vector<int>* T) {
    // update distribution for every net on BaseCell
    for (auto net_it = _cell_array[BaseCell->_ID].begin(); net_it != _cell_array[BaseCell->_ID].end(); ++net_it) {
        --(F->at(*net_it));
        ++(T->at(*net_it));
    }
}

bool Solver::update_max_gain_pointer() {
    for (int gain = _Pmax; gain >= -_Pmax; --gain) {
        if (_Bucket_A[gain].size()) {
            _maxGainPtr = &_Bucket_A[gain];
            _maxGainPtr->reset_pick_candidate();
            return true;
        }
        if (_Bucket_B[gain].size()) {
            _maxGainPtr = &_Bucket_B[gain];
            _maxGainPtr->reset_pick_candidate();
            return true;
        }
    }
    // cerr << "\tno more cells left in both buckets, process terminating." << endl;
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
    _gain_history.push_back(max_gain);
    if (max_gain > 0) {
        cerr << "\t> Iteration Gain: " << max_gain << endl;
        this->apply_change(k);
        _cutsize -= max_gain;
        return true;
    }
    else return false;
}

void Solver::apply_change(int k) {
    cerr << "\t> swapping cells" << endl;
    for (unsigned int i = k+1; i < _cell_gain_pairs.size(); ++i) {
        _cell_gain_pairs[i].first->change_group();
    } cerr << endl;
    this->initBucketList();
}

void Solver::assign_basecell(Cell* c, Group g) {
    BaseCell = c;
    _BaseCell_F = g;
}

void Solver::print_hisotry() const {
    cerr << endl << "[Gain History]" << endl;
    for (unsigned int i = 0; i < _gain_history.size(); ++i) {
        cerr << "\t> Iteration " << i+1 << ", Gain " << _gain_history[i] << endl;
    }
}

void Solver::dump(ostream& os) {
    cerr << "[Dumping result]" << endl;
    os << "Cutsize = " << _cutsize << endl;
    os << "G1 " << _Bucket_A.size() << endl;
    _Bucket_A.dump(os);
    os << "G2 " << _Bucket_B.size() << endl;
    _Bucket_B.dump(os);
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
