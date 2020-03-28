#ifndef DEFORMINFO_H
#define DEFORMINFO_H
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>


#define edl '\n'
#define FOR(i,l,r) for (int i=l; i<r; i++)
#define FORV(i,v) for (int i=0; i<(v).size(); i++)

using std::endl;

class CPHandle {
public :

    class CPG {
    public :
        std::vector<int> ids;
        std::vector<Eigen::Vector3d> cps;
        void clear() {ids.clear();cps.clear();}
    };

    std::vector<CPG> data;


    void Load(std::istream &cin) {
        data.clear();
        int n;
        cin >> n;
        CPG cpg;
        for (int i=0; i<n; i++) {
            cpg.clear();
            int m;
            cin >> m;
            for (int j=0; j<m; j++) {
                int id;
                cin >> id;
                cpg.ids.push_back(id);
            }
            for (int j=0; j<m; j++) {
                Eigen::Vector3d v;
                cin >> v(0) >> v(1) >> v(2);
                cpg.cps.push_back(v);
            }
            data.push_back(cpg);
        }
    }

    void Write(std::ostream &cout) {
        cout << data.size() << edl;
        FORV(i,data) {
            cout << data[i].ids.size() << edl;
            FORV(j,data[i].ids)
                    cout << data[i].ids[j] << edl;
            FORV(j,data[i].ids)
                    cout << data[i].cps[j].transpose() << edl;
            cout << endl;
        }
    }

    void Write(std::vector<double> &v) {
        FORV(i,data) {
            FORV(j,data[i].ids) {
                v.push_back(data[i].cps[j](0));
                v.push_back(data[i].cps[j](1));
                v.push_back(data[i].cps[j](2));
            }
        }
    }

    void Load(double *v) {
        FORV(i,data) {
            FORV(j,data[i].ids) {
                data[i].cps[j](0) = *(v++);
                data[i].cps[j](1) = *(v++);
                data[i].cps[j](2) = *(v++);
            }
        }
    }
};

template <class T>
void VectorPush(std::vector<T> &a, std::vector<T> &b) {
    for (int i=0; i<b.size(); i++) a.push_back(b[i]);
}

class DeformInfo {
public:
    std::vector<std::string> modelName;
    std::vector<double> weight;
    CPHandle handle;
    std::string prefix;

    static const char name1[];
    static const char name2[];

    void Load(std::string prefix) {
        this->prefix = prefix;
        std::ifstream cin((prefix + name1).c_str());
        handle.Load(std::ifstream((prefix + name2).c_str()));

        weight.clear();
        modelName.clear();

        int n;
        cin >> n;
        FOR(i,0,n) {
            double x;
            cin >> x;
            weight.push_back(x);
        }
        std::string x;
        std::getline(cin, x);
        FOR(i,0,n) {
            std::getline(cin, x);
            modelName.push_back(x);
        }
    }

    void Write() {
        std::ofstream cout((prefix + name1).c_str());
        handle.Write(std::ofstream((prefix + name2).c_str()));
        cout << modelName.size() << endl;
        FOR(i,0,modelName.size()) {
            cout << weight[i] << " ";
        }
        cout << endl;
        FOR(i,0,modelName.size()) {
            cout << modelName[i] << endl;
        }
    }

    void Write(std::vector<double> &v) {
        VectorPush(v, weight);
        handle.Write(v);
    }

    void Load(double *v) {
        for (int i=0; i<weight.size(); i++)
            weight[i] = *(v++);
        handle.Load(v);
    }
};

std::string ItoS(int i);

std::string operator +(std::string s, int i);
std::string operator +(int i, std::string s);
#endif // DEFORMINFO_H
