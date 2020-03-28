#include "DeformCaller.h"
#include "ARAPDeform.h"
#include <tuple>
#include <queue>

using namespace std;
using namespace Eigen;


void DVCaller::calc() {
    static bool isFirst = true;
    static ARAPDeform* deform;
    //static DMEngine eng;
    static DTriMesh result;

    if (isFirst) {
        isFirst = false;
        vector<DTriMesh*> v;
        for (int i=0; i<this->info.modelName.size(); i++) {
            DTriMesh *mesh = new DTriMesh();
            OpenMesh::IO::read_mesh(*mesh, info.modelName[i].c_str());
            v.push_back(mesh);
        }
        deform = new T2_ARAPDeform(*v[0], v);
        result = *v[0];
    }
    FeatureVector fv = deform->fvs[0];
    fv.isConst.resize(fv.s.size(), false);
    auto &d = info.handle.data;
    for (int i=0; i<d.size(); i++)
        for (int j=0; j<d[i].cps.size(); i++) {
            fv.isConst[d[i].ids[j]] = true;
            fv.constPoint[d[i].ids[j]] = d[i].cps[j];
        }
    string objname = info.prefix + ".obj";
    deform->solve2(fv, info.weight, result);
    OpenMesh::IO::write_mesh(result, objname.c_str());
}

void makeBarHandle() {
    DTriMesh mesh;
    OpenMesh::IO::read_mesh(mesh, "E:/SIGA2014/dataset/bar/bars/2.obj");
    double a = 9.358 - 0.001;
    double b = 0.415 + 0.001;

    CPHandle handle;
    handle.data.resize(2);

    for (int i=0; i<mesh.n_vertices(); i++) {
        Vector3d v = OtoE(mesh.point(DTriMesh::VertexHandle(i)));
        if (v(1) > a) {
            handle.data[0].ids.push_back(i);
            handle.data[0].cps.push_back(v);
        }
        if (v(1) < b) {
            handle.data[1].ids.push_back(i);
            handle.data[1].cps.push_back(v);
        }
    }
    cout << handle.data[0].ids.size() << endl;
    cout << handle.data[1].ids.size() << endl;
    handle.Write(ofstream("E:/SIGA2014/dataset/bar/bars/handle/c3.txt"));
}

#define FORE(i,v) for (decltype(v.begin()) i = v.begin(); i!=v.end(); i++)

void CallMainPre() {
    string objname = "E:\\SIGA2014\\dataset\\faces\\1.obj";
    map<int,bool> has;
    ifstream cin(objname);
    char x;
    string ha;
    getline(cin,ha);
    cout << ha;
    vector< tuple<int,int,int> > fs;
    ofstream cout("result.obj");
    int t=1;
    map<int,Vector3d> pp;
    while (cin>>x) {
        if (x != 'f') {
            double a,b,c;
            cin >> a >> b >> c;
            pp[t++] = Vector3d(a,b,c);
        } else {
            int a,b,c;
            cin >> a >> x >> x >> b >> x >> x >> c >> x >> x;
            has[a] = has[b] = has[c] = true;
            fs.push_back(make_tuple(a,b,c));
        }
    }
    map<int,int> rid;
    t=1;
    FORE(x,has) {
        int i = x->first;
        rid[i] = t++;
    }
    FORE(x,has) {
        int i = x->first;
        cout << "v " << pp[i].transpose() << endl;
    }
    FORE(x,fs) {
        cout << "f " << rid[get<0>(*x)] << " " << rid[get<1>(*x)] << " " << rid[get<2>(*x)] << endl;
    }
}

void CallMain() {
    string objname = "E:\\liangdun\\example\\dinos\\4.obj";
    cout << objname << endl;
    int root = 1;
    map<int,bool> has;
    ifstream cin(objname);
    string x;
    string ha;
    //getline(cin,ha);
    //cout << ha;
    vector< tuple<int,int,int,int> > fs;
    vector< tuple<int,int,int,int> > uvs;
    ofstream cout("hasuv.obj");
    int t=1;
    map<int,Vector3d> pp;
    map<int, map<int,bool> > e;
    while (cin>>x) {
        if (x == "v") {
            double a,b,c;
            cin >> a >> b >> c;
            pp[t++] = Vector3d(a,b,c);
        } else
        if (x == "f") {
            int a,b,c,d;
            int ua,ub,uc,ud;
            char tp;
            int tp2;
            cin >> a;
            cin >> tp >> ua >> tp >> tp2;
            cin >> b;
            cin >> tp >> ub >> tp >> tp2;
            cin >> c;
            cin >> tp >> uc >> tp >> tp2;
            cin >> d;
            cin >> tp >> ud >> tp >> tp2;
            e[a][b] = true;
            e[b][a] = true;
            e[a][c] = true;
            e[c][a] = true;
            e[b][c] = true;
            e[c][b] = true;
            e[c][d] = true;
            e[d][c] = true;
            fs.push_back(make_tuple(a,b,c,d));
            uvs.push_back(make_tuple(ua,ub,uc,ud));
        } else {
            getline(cin,x);
        }
    }
    queue<int> qu;
    qu.push(root);
    has[root] = true;
    while (!qu.empty()) {
        int i = qu.front();
        qu.pop();
        FORE(x,e[i]) {
            int j = x->first;
            if (!has.count(j)) {
                has[j] = true;
                qu.push(j);
            }
        }
    }
    map<int,int> rid;
    t=1;
    FORE(x,has) {
        int i = x->first;
        rid[i] = t++;
    }
    FORE(x,has) {
        int i = x->first;
        cout << "v " << pp[i].transpose() << endl;
    }
    int dd=-1;
    FORE(x,fs) {
        dd++;
        if (rid[get<0>(*x)] == 0) continue;
        cout << "f " <<
                rid[get<0>(*x)] << "/" << get<0>(uvs[dd]) << " " <<
                rid[get<1>(*x)] << "/" << get<1>(uvs[dd]) << " " <<
                rid[get<2>(*x)] << "/" << get<2>(uvs[dd]) << " " <<
                rid[get<3>(*x)] << "/" << get<3>(uvs[dd]) << " " <<
                endl;
    }
}
