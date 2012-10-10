/*  
    This file is part of HMMLab.

    HMMLab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HMMLab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HMMLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <gvc.h>
#include "vmlib.h"
#include "data_structures.h"
using namespace std;


#ifndef HMMLABLIB_H
#define HMMLABLIB_H

#define HTK_FORMAT "htk"
#define XML_FORMAT "xml"

class ModelSet;
class Model;
class State;
class Stream;
class Gaussian;
class SVector;
class SMatrix;
class HMMLab_Object;

string gettag(istream&);
void init();
string execute(string cmd);

enum hmmlab_types {
    MODELSET,
    MODEL,
    STATE,
    STREAM,
    GAUSSIAN,
    SVECTOR,
    TRANSMATRIX,
    SMATRIX,
    SHARED,
    HMMLAB_OBJECT
};

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

    virtual void __del__() {
        delete this;
    };
};

class HMMLab_Object : public RCObj
{
public:
    string name;
    hmmlab_types type;

    HMMLab_Object(): RCObj(), name(""), type(HMMLAB_OBJECT)  {};
    HMMLab_Object(string n, hmmlab_types t): RCObj(), name(n), type(t) {};
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
    void save(ostream&, const char*);
public:
    SVector(string, ModelSet*, int, double);

    friend class Gaussian;
    friend class ModelSet;
};


class SMatrix : public Shared, public Matrix
{
    void save(ostream&, const char*, bool);
public:
    SMatrix(string, ModelSet*, int, int, double);

    friend class Gaussian;
    friend class ModelSet;
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
    friend class ModelSet;
};


class Stream : public Shared
{
    void load(istream&, const char*, int, int);
    void save(ostream&, const char*);
public:
    int index_distribution;
    List<Gaussian*> gaussians;
    List<double> gaussians_weights;

    Stream(string, ModelSet*, int);
    Stream(string, ModelSet*, int, List<Gaussian*>, List<double>);
    ~Stream();

    void add_gaussian(Gaussian*, double);
    void remove_gaussian(Gaussian*);
    void remove_gaussian(int);
    double get_gaussian_weight(Gaussian*);

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
    List<List<Gaussian*>* > get_gaussians();

    friend class Model;
    friend class ModelSet;
};


class TransMatrix : public Shared
{
    void load(istream&, const char*, unsigned int);
    void save(ostream&, const char*, unsigned int);
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
    List<List<Gaussian*>* > get_gaussians();

    friend class ModelSet;
};


class ModelSet : public HMMLab_Object
{
    List<Vector* >* get_positions(graph_t*, unsigned int, const char*);
    List<Vector* >* translate_positions(List<Vector* >*);
    List<List<Vector* >* > orig_pos_data, data;
    List<List<double>* > edge_len;
    double screen_width, screen_height, orig_width, orig_heigth;

    void load(istream&, const char*);
    void save(ostream&, const char*);
    graph_t* layout_graph(unsigned int, GVC_t*);
    graph_t* layout_graph(unsigned int, GVC_t*, List<Vector*> gaussians_m);
    Vector* get_pos(graph_t*, char*);
    void add_data(List<Vector*>);

public:
    List<List<Vector*>*> pos_data;
    unsigned int dimension;
    unsigned int streams_size;
    List<unsigned int> streams_distribution;
    List<Model*> models;
    Dict<string, HMMLab_Object* > objects_dict;
    List<int> vecsize_tags;

    ModelSet();
    ModelSet(string);
    ModelSet(string, const char*);
    ~ModelSet();
    void destroy();
    void __del__();
    void save(const char*, const char*);

    void set_wh(double, double);
    List<List<Vector*>* > get_positions(List<List<Gaussian*>* >);
    List<List<Vector*>* > get_positions(double, double, List<List<Gaussian*>* >);
    void load_data(string);

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
