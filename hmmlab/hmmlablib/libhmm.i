%module(docstring = "Module for creating and manipulating HMMs") libhmm
%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%allowexception;

%{
#define SWIG_FILE_WITH_INIT
#include "libhmm.h"
using namespace std;
%}

%feature("ref")   RCObj "$this->inc_ref_num();"
%feature("unref") RCObj "$this->dec_ref_num();"

%include "data_structures.h"

%define vec_list_wrap(postfix, T)
    %template(vector ## postfix) std::vector<T >;
    %template(List ## postfix) List<T >;
%enddef

vec_list_wrap(Int, int)
vec_list_wrap(Double, double)
vec_list_wrap(Model, Model*)
vec_list_wrap(State, State*)
vec_list_wrap(Stream, Stream*)
vec_list_wrap(Gaussian, Gaussian*)
vec_list_wrap(Vector, Vector*)
vec_list_wrap(SVector, SVector*)
vec_list_wrap(ListDouble, List<double>*)
vec_list_wrap(ListListDouble, List<List<double>* >)
vec_list_wrap(HMMLab_Object, HMMLab_Object*)

%define map_dict_wrap(postfix, T, U)
    %template(map ## postfix) std::map<T , U >;
    %template(Dict ## postfix) Dict<T , U >;
%enddef

map_dict_wrap(StringHMMLab_Object, std::string, HMMLab_Object * )

%include "data_structures.h"

%rename(assign) *::operator=;
%rename(__getitem__) *::operator[];

%typemap(in)(ListInt)
{
    if(!PyList_Check($input)) {
        PyErr_SetString(PyExc_ValueError, "Expected a list");
        return NULL;
    }
    int len = PyList_Size($input);
    for(int i = 0; i < len; i++) {
        PyObject* o = PyList_GetItem($input, i);
        if(!PyLong_Check(o)) {
            PyErr_SetString(PyExc_ValueError, "Expected an integer");
            return NULL;
        }
        $1.append(PyLong_AsLong(o));
    }
}

%typemap(typecheck)(ListInt)
{
    $1 = PyList_Check($input) ? 1 : 0;
}

%typemap(in)(ListDouble)
{
    if(!PyList_Check($input)) {
        PyErr_SetString(PyExc_ValueError, "Expected a list");
        return NULL;
    }
    int len = PyList_Size($input);
    for(int i = 0; i < len; i++) {
        PyObject* o = PyList_GetItem($input, i);
        if(!PyFloat_Check(o)) {
            PyErr_SetString(PyExc_ValueError, "Expected a list float");
            return NULL;
        }
        $1.append(PyFloat_AsDouble(o));
    }
}

%typemap(typecheck)(ListDouble)
{
    $1 = PyList_Check($input) ? 1 : 0;
}

%init {
    init();
}

%include "data_structures.h"
%include "vmlib.h"

class ModelSet;

class RCObj
{
public:
    int ref_num;
    RCObj(){};
    virtual ~RCObj() {};

    virtual void inc_ref_num();
    virtual void dec_ref_num();
};

class HMMLab_Object : public RCObj
{
public:
    string name;

    HMMLab_Object(): name(""), type(HMMLAB_OBJECT) {};
    HMMLab_Object(string n, hmmlab_types t): name(n), type(t) {};
};


/* Trieda od ktorej budu vsetky zdielatelne dedit */
class Shared : public HMMLab_Object
{
public:
    ModelSet* modelset;

    Shared(string, hmmlab_types, ModelSet*);

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
public:
    int index_distribution;
    double gconst;
    SVector* mean;
    SMatrix* covariance;

    Gaussian(string, ModelSet*, int, double);
    Gaussian(string, ModelSet*, int, double, SVector*, SMatrix*);
    ~Gaussian();

};

class Stream : public Shared
{
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
public:
    int x, y;
    List<Stream*> streams;
    List<double> stream_weights;

    State(string, ModelSet*);
    State(string, ModelSet*, List<double>);
    State(string, ModelSet*, List<Stream*>, List<double>);
    ~State();
};

class TransMatrix : public Shared
{
    List<List<List<double> * >* > matrix;
public:
    TransMatrix(string, ModelSet*, int, double);
    ~TransMatrix();

    double operator()(unsigned int, int);
    void operator()(unsigned int, int, double);
    void remove(int);
    void remove_matrix(int);
    void add_matrix(TransMatrix&);
};


class Model : public HMMLab_Object
{
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
};


class ModelSet : public HMMLab_Object
{
public:
    int dimension;
    int streams_size;
    List<int> streams_distribution;
    List<Model*> models;
    Dict<string, HMMLab_Object* > objects_dict;
    List<int> vecsize_tags;

    ModelSet(string);
    ModelSet(string, const char*);
    ~ModelSet();
    void destroy();
    void save(const char*, const char*);

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
