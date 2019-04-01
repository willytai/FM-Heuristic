#ifndef SOLVER_H
#define SOLVER_H

#include "Dlist.h"
#include <vector>

using namespace std;

class Solver{
public:
    Solver() :  _balance_degree(0), _Pmax(0), _numCell(0), _numNet(0), _cutsize(0) {}
    ~Solver() {}

    // Top API
    void read(char* filename);
    void solve();
    void dump(ostream& os);

    void initPartition();
    void initBucketList();
    void construct_balance_criterion();
    void moveMaxGainCell();
    void move_max_gain_pointer(Group ref);
    void update_net_distribution(vector<int>* F, vector<int>* T);
    void apply_change(int k);
    void assign_basecell(Cell* c, Group g);
    void print_hisotry() const;
    bool update_gain();
    bool update_max_gain_pointer();
    bool compute_max_gain();
    bool balance_checking(Cell* c);

#ifdef DEBUG_MODE
    void debug_net_dist() const;
#endif

private:
    vector<vector<int> >    _cell_array;
    vector<vector<int> >    _net_array;
    vector<Cell*>           _cell_ptr;

    Bucket                  _Bucket_A;
    Bucket                  _Bucket_B;
    Dlist*                  _maxGainPtr;

    vector<int>             _NetADistribution;
    vector<int>             _NetBDistribution;
    
    Cell*                   BaseCell;
    Group                   _BaseCell_F;

    vector<pair<Cell*, int> > _cell_gain_pairs;

    double                  _balance_degree;

    double                  _max_limit;
    double                  _min_limit;
    
    int                     _Pmax;
    int                     _numCell;
    int                     _numNet;
    int                     _cutsize;
    
    // dubug
    vector<int>             _gain_history;
};

#endif /* SOLVER_H */
