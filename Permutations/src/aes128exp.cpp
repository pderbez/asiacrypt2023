#include <iostream>
#include <vector>
#include <set>
#include "CustomCallback.hpp"
#include "SysOfEqs.hpp"

//#include "/opt/gurobi/linux64/include/gurobi_c++.h"
//#include "/Library/gurobi952/macos_universal2/include/gurobi_c++.h"

using namespace std;

void printPerm(vector<int> & P) {
    for (unsigned i = 0; i < 16; ++i){
        if(P[i] == 16) cout << "* ";
        else cout << P[i] << " ";
    }
    cout << " " << endl;
}

vector<vector<unsigned>> modelAES128(int R, vector<int> & KPerm, int nrSboxesWanted, GRBEnv & env) {

    GRBModel model = GRBModel(env);

    auto const & size_perm = KPerm.size();
    auto const & nb_cols_perm = size_perm/4;


    // X = 3r + 0
    // S(X) = -(3*r + 0)
    // Y (after MC) = 3*r + 1
    // K = 3*r + 2

    // definition of variables
    vector<GRBVar> dX ((3*R + 1)*16); // before ARK

    for (unsigned i = 0; i < (3*R+1)*16; ++i) {
        dX[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }

    for (unsigned i = 0; i < 3*16; ++i) {
      model.addConstr(dX[i] == 0);
    }


    // Contraints for permutation-based key schedule
    // number of active bytes per key must be the same
    // {
    //   GRBLinExpr e = 0;
    //   for (unsigned i = 0; i < 16; ++i) e += dX[16*(3*1 + 2) + i];
    //   for (unsigned r = 2; r < R; ++r) {
    //     GRBLinExpr e1 = 0;
    //     for (unsigned i = 0; i < 16; ++i) e1 += dX[16*(3*r + 2) + i];
    //     model.addConstr(e1 == e);
    //   }
    // }

    vector<vector<unsigned>> subkeys;

    {
      vector<unsigned> v (size_perm);
      for (unsigned x = 0; x < size_perm; ++x) v[x] = x;
      subkeys.emplace_back(v);
      while (size_perm*subkeys.size() < 16*R) {
        vector<unsigned> vv (size_perm);
        for (unsigned i = 0; i < size_perm; ++i) vv[KPerm[i]] = v[i];
        subkeys.emplace_back(vv);
        v = move(vv);
      }
    }

    for (unsigned c = 4; c < 4*R; ++c) {
      unsigned r_roundk = c/4;
      unsigned c_roundk = c%4;
      unsigned r_subk = (c-4)/nb_cols_perm;
      unsigned c_subk = (c-4)%nb_cols_perm;
      if (r_subk == 0) continue;
      for (unsigned l = 0; l < 4; ++l) {
        unsigned x = subkeys[r_subk][4*c_subk + l];
        unsigned rx = 1;
        unsigned cx = x%nb_cols_perm;
        unsigned lx = x/nb_cols_perm;
        while (cx >= 4) {cx -= 4; rx += 1;}
        if (rx < R) model.addConstr(dX[16*(3*rx + 2) + 4*lx + cx] == dX[16*(3*r_roundk + 2) + 4*l + c_roundk]);
      }
    }


    for (unsigned r = 1; r < R; ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        GRBLinExpr e = 0;
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 1)  + c  + 4*i];
        for (unsigned i = 0; i < 4; ++i) e += dX[16*(3*r + 0) + ((c + i)%4) + 4*i];
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        model.addConstr(e <= 8*f);
        model.addConstr(e >= 5*f);
      }
    }

    for (unsigned r = 1; r < R ; ++r) {
      //ARK (sX[r+1], kX[r+1] and zX[r+1])
      for (unsigned i = 0; i < 16; ++i) {
        model.addConstr(1-dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + 1-dX[16*(3*r + 2)  + i] + dX[16*(3*(r+1))  + i] >= 1);
        model.addConstr(dX[16*(3*r + 1)  + i] + dX[16*(3*r + 2)  + i] + 1-dX[16*(3*(r+1))  + i] >= 1);
      }
    }

    GRBLinExpr obj = 0;
    for (unsigned r = 1; r <= R; ++r) {
      for (unsigned i = 0; i < 16; ++i) obj += dX[16*(3*r) + i];
    }

    GRBLinExpr obj2 = 0;
    //for (unsigned i = 0; i < 16; ++i) obj2 += dX[16*(3*1 + 2) + i];
    for (unsigned i = 0; i < size_perm; ++i) {
      unsigned r = 1;
      unsigned c = i%nb_cols_perm;
      unsigned l = i/nb_cols_perm;
      while (c >= 4) {c -= 4; r += 1;}
      if (r < R) obj2 += dX[16*(3*r + 2) + 4*l + c];
    }

    model.setObjective(obj2, GRB_MINIMIZE);

    //model.addConstr(obj2 == 2);




    /*GRBLinExpr obj3 = 0;

    for (unsigned r = 1; r < R; ++r) {
      for (unsigned c = 0; c < 4; ++c) {
        GRBLinExpr e1 = 0;
        GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        GRBVar g = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        for (unsigned i = 0; i < 4; ++i) e1 += dX[16*(3*r + 0) + ((c + i)%4) + 4*i];
        model.addConstr(e1 <= 4*f + 4*g);
        GRBLinExpr e2 = 0;
        for (unsigned i = 0; i < 4; ++i) e2 += dX[16*(3*(r+1) + 0)  + c  + 4*i];
        model.addConstr(e1 + e2 >= 4*(1-f) - 4*(1-g));
        obj3 += f;
      }
    }

    model.setObjective(obj3, GRB_MINIMIZE);*/

    model.addConstr(obj >= 1);
    model.addConstr(obj <= nrSboxesWanted-1);

    auto mat = AES128eqs(R, KPerm.data(), subkeys);

    mycallback cb ((3*R+1)*16, dX.data(), mat);
    model.setCallback(&cb);
    //model.set(GRB_IntParam_OutputFlag , (R == 6) ? 1 : 0);
    model.set(GRB_IntParam_LazyConstraints , 1);
    // model.set(GRB_IntParam_PoolSolutions, 2000000);
    // model.set(GRB_DoubleParam_PoolGap, 0.001);
    // model.set(GRB_IntParam_PoolSearchMode, 2);

    //printPerm(KPerm);

    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);


    vector<vector<unsigned>> v_res;
    if (nSolutions > 0) {
      for (auto e = 0; e < 1; e++) {
          model.set(GRB_IntParam_SolutionNumber, e);

          vector<unsigned> res;
          unsigned c_sub = 0;
          for (unsigned r = 1; r < R; ++r) {
            for (unsigned c = 0; c < 4; ++c) {
              for (unsigned l = 0; l < 4; ++l) {
                if (dX[16*(3*r + 2) + 4*l + c].get(GRB_DoubleAttr_Xn) > 0.5) res.emplace_back(4*l + c_sub);
              }
              c_sub += 1;
              if (c_sub == nb_cols_perm) {
                c_sub = 0;
                v_res.emplace_back(res);
                res.clear();
              }
            }
          }




      }
    }


    return v_res;

}


void searchMILP(vector<pair<unsigned, unsigned>> const & myconstraints, int size_perm) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  vector<int> KPerm(size_perm, size_perm);

  GRBModel model = GRBModel(env);
  model.set(GRB_IntParam_OutputFlag , 0);

  unsigned max_R = 0;
  for (auto const & p : myconstraints) {
    if (max_R < p.first) max_R = p.first;
  }

  if (max_R > 2) max_R -= 1;

  vector<vector<vector<GRBVar>>> P (max_R, vector<vector<GRBVar>> (size_perm, vector<GRBVar> (size_perm)));

  for (unsigned r = 1; r < max_R; ++r) {
    for (unsigned i = 0; i < size_perm; ++i) {
      for (unsigned j = 0; j < size_perm; ++j) P[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
    for (unsigned i = 0; i < size_perm; ++i) {
      GRBLinExpr e = 0;
      for (unsigned j = 0; j < size_perm; ++j) e += P[r][i][j];
      model.addConstr(e == 1);
    }

    for (unsigned i = 0; i < size_perm; ++i) {
      GRBLinExpr e = 0;
      for (unsigned j = 0; j < size_perm; ++j) e += P[r][j][i];
      model.addConstr(e == 1);
    }
  }


  auto const n_cols_perm = size_perm/4;

  // remove symetries
  for (unsigned j0 = 0; j0 < size_perm; ++j0) {
    for (unsigned i = 1; i < 4; ++i) {
      for (unsigned j = 0; j < size_perm; ++j) {
        unsigned jj = n_cols_perm*(j/n_cols_perm) + (j%n_cols_perm - i + n_cols_perm)%n_cols_perm;
        if (jj < j0) model.addConstr(P[1][i][j] <= 1 - P[1][0][j0]);
      }
    }
  }

  for (unsigned r = 2; r < (4*max_R+3)/n_cols_perm; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      for (unsigned j = 0; j < 16; ++j) {
        for (unsigned k = 0; k < 16; ++k) model.addConstr(P[r][i][k] + 1 >= P[r-1][i][j] + P[1][j][k]);
      }
    }
  }




  bool found = false;
  while (!found) {
    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);
    if (nSolutions == 0) {cout << "pas de solutions" << endl; getchar();}
    model.set(GRB_IntParam_SolutionNumber, 0);
    for (unsigned i = 0; i < size_perm; ++i) {
      unsigned j = 0;
      while (P[1][i][j].get(GRB_DoubleAttr_Xn) < 0.5) ++j;
      KPerm[i] = j;
    }
    for (unsigned i = 0; i < size_perm; ++i) cout << KPerm[i] << " ";
    cout << endl;

    unsigned j = 0;

    unsigned R = myconstraints[0].first;
    unsigned nbsboxes = myconstraints[0].second;

    auto v_res = modelAES128(R, KPerm, nbsboxes, env);
    while (v_res.empty() && ++j < myconstraints.size()) {
    	R = myconstraints[j].first;
    	nbsboxes = myconstraints[j].second;
    	v_res = modelAES128(R, KPerm, nbsboxes, env);
    }
    if (v_res.empty()) {
      cout << "yeah!!" << endl;
      for (unsigned i = 0; i < size_perm; ++i) cout << KPerm[i] << " ";
      getchar();
      GRBLinExpr e = 0;
      for (unsigned i = 0; i < size_perm; ++i) e += P[1][i][KPerm[i]];
      model.addConstr(e <= size_perm-1);
    }
    else {
      cout << "v_res: " << v_res.size() << endl;
      unsigned bound = 0;
      GRBLinExpr e = 0;
      for (unsigned r = 1; r < v_res.size(); ++r) {
        for (auto i1 : v_res[r-1]) {
          for (auto i2 : v_res[r]) e += P[1][i1][i2];
          bound += 1;
        }
      }
      model.addConstr(e <= bound - 1);

    }
  }
}





int main(int argc, char const *argv[]) {
    searchMILP(vector<pair<unsigned, unsigned>>({make_pair(8, 22)}), 32);






    return 0;
}
