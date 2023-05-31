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
    {
      GRBLinExpr e = 0;
      for (unsigned i = 0; i < 16; ++i) e += dX[16*(3*1 + 2) + i];
      for (unsigned r = 2; r < R; ++r) {
        GRBLinExpr e1 = 0;
        for (unsigned i = 0; i < 16; ++i) e1 += dX[16*(3*r + 2) + i];
        model.addConstr(e1 == e);
      }
    }

    for (unsigned i = 0; i < 16; ++i) {
      unsigned x = i;
      for (unsigned r = 2; r < R; ++r) {
          model.addConstr(dX[16*(3*r + 2) + KPerm[x]] == dX[16*(3*(r-1) + 2) + x]);
          x = KPerm[x];
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
    for (unsigned i = 0; i < 16; ++i) obj2 += dX[16*(3*1 + 2) + i];

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

    auto mat = AES128eqs(R, KPerm.data());

    mycallback cb ((3*R+1)*16, dX.data(), mat);
    model.setCallback(&cb);
    model.set(GRB_IntParam_OutputFlag , (R == 6) ? 1 : 0);
    model.set(GRB_IntParam_LazyConstraints , 1);
    model.set(GRB_IntParam_PoolSolutions, 2000000);
    model.set(GRB_DoubleParam_PoolGap, 0.001);
    model.set(GRB_IntParam_PoolSearchMode, 2);

    //printPerm(KPerm);

    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);


    vector<vector<unsigned>> v_res;

    for (auto e = 0; e < nSolutions; e++) {
        model.set(GRB_IntParam_SolutionNumber, e);

        vector<unsigned> res;
        for (unsigned i = 0; i < 16; ++i) {
            if (dX[16*(3*1 + 2) + i].get(GRB_DoubleAttr_Xn) > 0.5) res.emplace_back(i);
        }
        v_res.emplace_back(move(res));
    }

    return v_res;

}


void searchMILP(vector<pair<unsigned, unsigned>> const & myconstraints) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  vector<int> KPerm(16, 16);

  GRBModel model = GRBModel(env);
  model.set(GRB_IntParam_OutputFlag , 0);
  
  unsigned max_R = 0;
  for (auto const & p : myconstraints) {
    if (max_R < p.first) max_R = p.first;
  }
  
  if (max_R > 2) max_R -= 1;

  vector<vector<vector<GRBVar>>> P (max_R, vector<vector<GRBVar>> (16, vector<GRBVar> (16)));

  for (unsigned r = 1; r < max_R; ++r) {
    for (unsigned i = 0; i < 16; ++i) {
      for (unsigned j = 0; j < 16; ++j) P[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
    for (unsigned i = 0; i < 16; ++i) {
      GRBLinExpr e = 0;
      for (unsigned j = 0; j < 16; ++j) e += P[r][i][j];
      model.addConstr(e == 1);
    }

    for (unsigned i = 0; i < 16; ++i) {
      GRBLinExpr e = 0;
      for (unsigned j = 0; j < 16; ++j) e += P[r][j][i];
      model.addConstr(e == 1);
    }
  }


   for (unsigned c = 0; c < 4; ++c) {
     for (unsigned l1 = 0; l1 < 4; ++l1) {
       for (unsigned l2 = 0; l2 < 4; ++l2) {
         model.addConstr(P[1][4*l1 + ((c+l1)%4)][4*l2 + c] == 0);
         //model.addConstr(P[1][4*l2 + c][4*l1 + ((c+l1)%4)] == 0);
       }
     }
   }

 


  for (unsigned j0 = 0; j0 < 16; ++j0) {
    for (unsigned i = 1; i < 4; ++i) {
      for (unsigned j = 0; j < 16; ++j) {
        unsigned jj = 4*(j/4) + (j%4 - i + 4)%4;
        if (jj < j0) model.addConstr(P[1][i][jj] <= 1 - P[1][0][j0]);
      }
    }
  }

  for (unsigned r = 2; r < max_R; ++r) {
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
    for (unsigned i = 0; i < 16; ++i) {
      unsigned j = 0;
      while (P[1][i][j].get(GRB_DoubleAttr_Xn) < 0.5) ++j;
      KPerm[i] = j;
    }
    for (unsigned i = 0; i < 16; ++i) cout << KPerm[i] << " ";
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
      for (unsigned i = 0; i < 16; ++i) cout << KPerm[i] << " ";
      getchar();
      GRBLinExpr e = 0;
      for (unsigned i = 0; i < 16; ++i) e += P[1][i][KPerm[i]];
      model.addConstr(e <= 15);
    }
    else {
      cout << "v_res: " << v_res.size() << endl;
      for (auto & res : v_res) {
        vector<vector<unsigned>> myvec;
        myvec.emplace_back(res);
        for (unsigned r = 1; r < R-1; ++r) {
          for (auto & x : res) x = KPerm[x];
          myvec.emplace_back(res);
        }
        vector<vector<unsigned>> myvec_col = myvec;
        for (auto & v : myvec_col) {
          unsigned cpt[4] = {0,0,0,0};
          for (auto x : v) cpt[x%4] += 1;
          for (auto & x : v) {
            if (cpt[x%4] > 1) x = 1;
            else x = 0;
          }
        }

        unsigned bound = 0;
        GRBLinExpr e = 0;
        for (unsigned r = 0; r < R-2; ++r) {
          for (auto x : myvec[r]) {
            for (auto y : myvec[r+1]) e += P[1][x][y];
            bound += 1;
          }
        }
        for (unsigned i = 0; i < res.size(); ++i) {
          for (unsigned r1 = 0; r1 < R-2; ++r1) {
            if (myvec_col[r1][i] == 0) continue;
            for (unsigned r2 = r1+1; r2 < R-1; ++r2) {
              if (myvec_col[r2][i] != 0) {
                e += P[r2-r1][myvec[r1][i]][myvec[r2][i]];
                bound += 1;
              }
            }
          }
        }
        model.addConstr(e <= bound - 1);
      }
    }
  }
}



bool modelAES128(int R, vector<int> & KPerm, int nrSboxesWanted, int active, GRBEnv & env) {

    GRBModel model = GRBModel(env);

    //model.set(GRB_IntParam_PoolSolutions, 2000000);
    //model.set(GRB_DoubleParam_PoolGap, 0.001);
    //model.set(GRB_IntParam_PoolSearchMode, 2);

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
    {
      GRBLinExpr e = 0;
      for (unsigned i = 0; i < 16; ++i) e += dX[16*(3*1 + 2) + i];
      for (unsigned r = 2; r < R; ++r) {
        GRBLinExpr e1 = 0;
        for (unsigned i = 0; i < 16; ++i) e1 += dX[16*(3*r + 2) + i];
        model.addConstr(e1 == e);
      }
    }

    for (unsigned i = 0; i < 16; ++i) {
      unsigned x = i;
      unsigned r = 1;
      while (r < R && x != 16) {++r; x = KPerm[x];}
      if (r == R) { // OK
        x = i;
        for (r = 2; r < R; ++r) {
            model.addConstr(dX[16*(3*r + 2) + KPerm[x]] == dX[16*(3*(r-1) + 2) + x]);
            x = KPerm[x];
        }
      }
      else { // NOK
        model.addConstr(dX[16*(3*1 + 2) + i] == 0);
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

    // path should involve the new guess
    {
      GRBLinExpr e = 0;
      for (unsigned r = 1; r < R; ++r) e += dX[16*(3*r + 2) + active];
      model.addConstr(e >= 1);
    }


    model.addConstr(obj >= 1);
    model.addConstr(obj <= nrSboxesWanted-1);

    auto mat = AES128eqs(R, KPerm.data());

    mycallback cb ((3*R+1)*16, dX.data(), mat);
    model.setCallback(&cb);
    model.set(GRB_IntParam_OutputFlag , 0);
    model.set(GRB_IntParam_LazyConstraints , 1);

    //printPerm(KPerm);

    model.optimize();

    auto nSolutions = model.get(GRB_IntAttr_SolCount);

    return (nSolutions != 0);
}

int PermComplete(unsigned *P){
    for(int i = 0; i<16; i++){
        if(P[i] == 16) return 0;
    }
    return 1;
}

unsigned firstAvailableElm(unsigned *P){
  /*  unsigned min = 16;
    for(int i = 0; i<16; i++){
        if((P[i] != 16) && (P[i] < min)) min = P[i];
    }
    return min;*/
}

void searchRec(unsigned R, int nrSboxesWanted, vector<int> & KPerm, vector<int> & images, unsigned n_images, unsigned n_forbidden, int x, unsigned length, unsigned min_length, GRBEnv & env) {

   if (KPerm[x] != 16) {
       if (length < min_length || (n_images > 0 && n_images < length)) return;
       // toujours évaluer le cycle
       if (modelAES128(R, KPerm, nrSboxesWanted, x, env)) return;

       // regarder si l'on peut commencer un autre cycle ou si on a tout fixé
       if (n_images == 0) {
         cout << "une permutation choisie: " << endl;
         for (const auto & y : KPerm) cout << y << " ";
         cout << endl;
         getchar();
       }
       else {
         if (n_forbidden == 0) {
           searchRec(R, nrSboxesWanted, KPerm, images, n_images, 0, images[0], 0, length, env);
         }
         else if (n_images >= length + 1) searchRec(R, nrSboxesWanted, KPerm, images, n_images, 0, images[0], 0, length+1, env);

         if (n_images >= 2*length + 1) {
           unsigned n_f = 1;
           if (n_images >= 2*length + 2) {
             while (n_f < n_forbidden && n_images - n_f >= length+1) {
               searchRec(R, nrSboxesWanted, KPerm, images, n_images, n_f, images[n_f], 0, length + 1, env);
               ++n_f;
             }
           }
           n_f = n_forbidden;
           while (n_images - n_f >= length) {
             searchRec(R, nrSboxesWanted, KPerm, images, n_images, n_f, images[n_f], 0, length, env);
             ++n_f;
           }
         }

       }
   }
   else {
       // regarder si l'on peut évaluer le début du cycle (si length est suffisamment grand)
       if ((length >= R-1) && modelAES128(R, KPerm, nrSboxesWanted, x, env)) return;
       // continuer le début du cycle
       for (unsigned i = n_forbidden; i < n_images; ++i) {
         KPerm[x] = images[i];
         images[i] = images[n_images - 1];
         searchRec(R, nrSboxesWanted, KPerm, images, n_images-1, n_forbidden, KPerm[x], length+1, min_length, env);
         images[i] = KPerm[x];
       }
       KPerm[x] = 16;
   }
}

void search(unsigned R, int nrSboxesWanted) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  vector<int> KPerm(16, 16);
  vector<int> images(16);
  for (int i = 0; i < 16; ++i) images[i] = i;
  searchRec(R, nrSboxesWanted, KPerm, images, 16, 0, 0, 0, 0, env);
  searchRec(R, nrSboxesWanted, KPerm, images, 16, 4, 4, 0, 0, env);
  searchRec(R, nrSboxesWanted, KPerm, images, 16, 8, 8, 0, 0, env);
  searchRec(R, nrSboxesWanted, KPerm, images, 16, 12, 12, 0, 0, env);
}

void searchCyclesRec(unsigned R, int nrSboxesWanted, vector<int> & KPerm, vector<int> & images, unsigned n_images, int x, unsigned length, unsigned length_obj, GRBEnv & env, bool & found) {
   if (KPerm[x] != 16) {
       if (length != length_obj) return;
       // toujours évaluer le cycle
       if (modelAES128(R, KPerm, nrSboxesWanted, x, env)) return ;

       int y = KPerm[x];
       cout << "\r" << x << " --> ";
       while (y != x) {
         cout << y << " --> ";
         y = KPerm[y];
       }
       cout << x << endl;

       found = true;

       if (length == 16) {cout << "permutation!!" << endl; getchar();}
   }
   else {
      if (length == length_obj) return;
       // regarder si l'on peut évaluer le début du cycle (si length est suffisamment grand)
       if ((length >= R-1) && modelAES128(R, KPerm, nrSboxesWanted, x, env)) return;
       // continuer le début du cycle
       for (unsigned i = 0; i < n_images; ++i) {
         KPerm[x] = images[i];
         images[i] = images[n_images - 1];
         searchCyclesRec(R, nrSboxesWanted, KPerm, images, n_images-1, KPerm[x], length+1, length_obj, env, found);
         images[i] = KPerm[x];
       }
       KPerm[x] = 16;
   }
}

void searchCycles(unsigned R, int nrSboxesWanted) {
  GRBEnv env = GRBEnv(true);
  env.set("LogFile", "mip1.log");
  env.start();
  vector<int> KPerm(16, 16);
  vector<int> images(16);
  for (int i = 0; i < 16; ++i) images[i] = 15-i;
  vector<bool> found (17, false);
  vector<bool> possible (17, false);
  possible[0] = true;
  for (unsigned l = 1; l <= 16; ++l) {
    bool found_cycle = false;
    if (l <= 8 || possible[16-l]) {
      searchCyclesRec(R, nrSboxesWanted, KPerm, images, 16, 0, 0, l, env, found_cycle);
      if (l <= 12) searchCyclesRec(R, nrSboxesWanted, KPerm, images, 12, 4, 0, l, env, found_cycle);
      if (l <= 8) searchCyclesRec(R, nrSboxesWanted, KPerm, images, 8, 8, 0, l, env, found_cycle);
      if (l <= 4) searchCyclesRec(R, nrSboxesWanted, KPerm, images, 4, 12, 0, l, env, found_cycle);
    }
    found[l] = found_cycle;
    if (found_cycle) {
      for (unsigned i = l; i <= 16; i += l) {
        for (unsigned j = 0; j <= 16; j += 1) {
          if (possible[j] && j + i <= 16) possible[j+i] = true;
        }
      }
    }
    cout << "cycles " << l << ": " << (found_cycle ? "possible" : "impossible") << endl;
  }

}



int main(int argc, char const *argv[]) {
    searchMILP({5, 15});
    //search(5, 17);
    //searchCycles(5, 18);





    return 0;
}
