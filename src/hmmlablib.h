#include <stdio.h>
#include <string.h>
#include <gvc.h>
#include <exception>
#include <fstream>
#include <iostream>
#include "vmlib.h"
#include "data_structures.h"
using namespace std;


#ifndef HMMLABLIB_H
#define HMMLABLIB_H

class ModelSet;
class Model;
class State;
class Stream;
class Gaussian;
class SVector;
class SMatrix;
class HMMLab_Object;

string gettag(ifstream&);
void init();

class RCObj
{
public:
    int ref_num;
    RCObj() : ref_num(0) {};
    virtual ~RCObj() {};

    virtual void inc_ref_num() {
        ref_num++;
    };
    virtual void dec_ref_num() {
        ref_num--;
        if(ref_num <= 0) {
            delete this;
        }
    };
};

class HMMLab_Object : public RCObj
{
public:
    string name;

    HMMLab_Object(): RCObj(), name("") {};
    HMMLab_Object(string n): RCObj(), name(n) {};
};


/* Trieda od ktorej budu vsetky zdielatelne dedit */
class Shared : public HMMLab_Object
{
public:
    ModelSet* modelset;

    Shared(string, ModelSet*);
    void inc_ref_num();
    void dec_ref_num();
};


class SVector : public Shared, public Vector
{
public:
    SVector(string, ModelSet*, int, double);
};


class SMatrix : public Shared, public Matrix
{
public:
    SMatrix(string, ModelSet*, int, int, double);
};


class Gaussian : public Shared
{
    void load(istream&, const char*);
    void save(ostream&, const char*);
public:
    int index_distribution;
    double gconst;
    SVector* mean;
    SMatrix* covariance;

    Gaussian(string, ModelSet*, int, double);
    Gaussian(string, ModelSet*, int, double, SVector*, SMatrix*);
    ~Gaussian();

    friend class Stream;
};


class Stream : public Shared
{
    void load(istream&, const char*, int, int);
    void save(ostream&, const char*);
    graph_t* layout_graph(GVC_t*);
    Vector* get_pos(graph_t*, char*);
    List<Vector* > get_positions(graph_t*, int, const char*);
public:
    int index_distribution;
    double screen_width, screen_height;
    List<Gaussian*> gaussians;
    List<double> gaussians_weights;
    List<Vector* > pos_data, pos_gauss, data;
    List<double> edge_len;

    Stream(string, ModelSet*, int);
    Stream(string, ModelSet*, int, List<Gaussian*>, List<double>);
    ~Stream();

    List<SVector*> gaussians_means();
    void add_gaussian(Gaussian*, double);
    void remove_gaussian(Gaussian*);
    void remove_gaussian(int);
    double get_gaussian_weight(Gaussian*);
    void add_data(List<Vector*>);
    void set_wh(double, double);

    friend class State;
};

class State : public Shared
{
    void load(istream&, const char*);
    void save(ostream&, const char*);
public:
    int x, y;
    List<Stream*> streams;
    List<double> stream_weights;

    State(string, ModelSet*);
    State(string, ModelSet*, List<double>);
    State(string, ModelSet*, List<Stream*>, List<double>);
    ~State();

    friend class Model;
    friend class ModelSet;
};


class TransMatrix : public Shared
{
    void load(istream&, const char*, int);
    void save(ostream&, const char*);
    List<List<List<double> * >* > matrix;
public:
    TransMatrix(string, ModelSet*, int, double);
    ~TransMatrix();

    double operator()(unsigned int, int);
    void operator()(unsigned int, int, double);
    void operator++();
    void remove(int);
    void remove_matrix(int);
    void add_matrix(TransMatrix&);

    friend class Model;
    friend class ModelSet;
};


class Model : public HMMLab_Object
{
    void load(istream&, const char*);
    void save(ostream&, const char*);
public:
    ModelSet* modelset;
    List<State*> states;
    TransMatrix* trans_mat;

    Model(string, ModelSet*);
    Model(string, ModelSet*, List<State*>, TransMatrix*);
    ~Model();

    void add_state(State*);
    void remove_state(State*);
    void remove_state(int);

    friend class ModelSet;
};


class ModelSet : public HMMLab_Object
{
    void load(istream&, const char*);
    void save(ostream&, const char*);
public:
    int dimension;
    int streams_size;
    List<int> streams_distribution;
    List<Model*> models;
    Dict<string, HMMLab_Object* > objects_dict;
    List<int> vecsize_tags;

    ModelSet();
    ModelSet(string);
    ModelSet(string, const char*);
    ~ModelSet();
    void destroy();

    void add_model(Model*);
    void remove_model(Model*);
    void remove_model(int);

    Model* get_model(string);
    State* get_state(string);
    Stream* get_stream(string);
    Gaussian* get_gaussian(string);
    SVector* get_svector(string);
    SMatrix* get_smatrix(string);
};

#endif
