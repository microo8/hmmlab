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
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <set>
#include <gvc.h>
#include "vmlib.h"
#include "data_structures.h"
#include "gnuplot_pipes.h"
using namespace std;


#ifndef HMMLABLIB_H
#define HMMLABLIB_H

#define HTK_FORMAT "htk"
#define XML_FORMAT "xml"
#define BORDER 10

class ModelSet;
class Model;
class State;
class Stream;
class Gaussian;
class SVector;
class SMatrix;
class HMMLab_Object;
class StreamArea;

string gettag(istream&);
void init();
string execute(string cmd);

struct point_len {
    unsigned int i;
    unsigned int j;
    double len;
    point_len(unsigned int ii, unsigned int jj, double l): i(ii), j(jj), len(l) {};
};

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
    void calc_gconst();
public:
    int index_distribution;
    double gconst;
    SVector* mean;
    SMatrix* covariance, *inv_covariance;
    List<unsigned int> my_data;

    Gaussian(string, ModelSet*, int, double);
    Gaussian(string, ModelSet*, int, double, SVector*, SMatrix*);
    ~Gaussian();

    double probability(Vector*);
    void divide();

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
    void select_gaussians(bool);
    void unselect_gaussians(bool);
    Gaussian* get_gaussian(unsigned int, bool);
    bool has_gaussian(Gaussian*);

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
    void select_gaussians();
    void unselect_gaussians();
    string create_image();

    friend class ModelSet;
};

class StreamArea
{
    double screen_width, screen_height; //velkost DrawArea
    double graph_width, graph_height; //velkost grafu
    double graph_prob_width, graph_prob_height; //velkost grafu s pravdepodobnostami
    double pca_width, pca_height; //sirka a vyska PCA dat
    double edge_len_multiplier; //prisposoby dlzky hran medzi datami, aby najmensia dlzka bola 1
    List<Vector*> data; //data pripadajuce na tento stream
    List<Vector*> last_pos_data; //pozicie na grafe v poslednom layoute
    List<Vector*> last_pos_data_pca; //pozicie dat v poslednom vypocitanom PCA
    List<Vector*> last_gauss_pos; //pozicie stredov gaussianov v poslednom
    List<Vector*> last_gauss_pos_pca; //pozicie stredov gaussianov v poslednom vypocitanom PCA
    List<Vector*> last_gauss_var_pca; //variancie gaussianov v PCA

    List<Vector*> last_pos_data_prob; //pozicie dat v poslednom layoute s dlzkami hran podla pravdepodobnosti
    List<Vector*> last_gauss_pos_prob; //poslednom stredov gaussianov v poslednom layoute s dlzkami hran podla pravdepodobnosti
    List<double> edge_len; //vzdialenosti medzi datami

    List<Vector* >* get_positions(graph_t*, unsigned int, const char*, bool prob);
    List<Vector* > translate_positions(List<Vector* >*);
    List<Vector* > translate_positions_prob(List<Vector* >*);
    List<Vector* > translate_pca_positions(List<Vector* >*, bool);
    graph_t* layout_graph(GVC_t*, bool);
    graph_t* layout_graph(GVC_t*, List<Vector*>);
    graph_t* layout_graph_prob(GVC_t*);
    Vector* get_pos(graph_t*, char*);

    void save_data_pos_2D(unsigned int, string);
    void save_data_pos_3D(unsigned int, unsigned int, string);

public:
    ModelSet* modelset;
    List<Vector*> pos_data; //pozicie translatovane na velkost DrawArea
    List<Vector*> pos_data_pca; //pozicie dat 2D PCA
    List<Vector*> pos_gaussians;
    List<Vector*> pos_gaussians_pca; //pozicie stredov gaussianov 2D PCA
    List<Vector*> pos_gaussians_var_pca; //variancie gaussianov v 2D PCA

    List<Vector*> pos_data_prob;
    List<Vector*> pos_gaussians_prob;

    set<Gaussian*> selected_gaussians;

    StreamArea(ModelSet*);
    ~StreamArea();
    Vector* get_data(unsigned int);
    void add_data(List<Vector*>);
    void set_wh(double, double);
    void reset_pos_gauss();
    void calc_pca();
    List<Vector*> get_data_2D(unsigned int, unsigned int);
    void calc_data_gauss();

    friend class ModelSet;
};


class ModelSet : public HMMLab_Object
{
    void load(istream&, const char*);
    void save(ostream&, const char*);
    void create_cfg(string);
    void add_data(List<Vector*>);
public:
    unsigned int dimension;
    unsigned int streams_size;
    List<unsigned int> streams_distribution;
    List<Model*> models;
    Dict<string, HMMLab_Object* > objects_dict;
    List<int> vecsize_tags;
    List<StreamArea*> stream_areas;

    ModelSet();
    ModelSet(string);
    ModelSet(string, const char*);
    ~ModelSet();
    void destroy();
    void save(const char*, const char*);

    void load_data(unsigned int, string*);
    void reset_pos_gauss();

    void add_model(Model*);
    void remove_model(Model*);
    void remove_model(int);

    bool is_selected(Model*, int);
    bool is_selected(Gaussian*);
    unsigned int selected_gaussians_count();
    unsigned int loaded_data_count();
    List<Model*> get_models_with_gaussian(Gaussian*);
    string get_unique_name(string);

    bool gauss_cluster(List<Gaussian*>, List<Vector*>);

    Model* get_model(string);
    State* get_state(string);
    Stream* get_stream(string);
    Gaussian* get_gaussian(string);
    SVector* get_svector(string);
    SMatrix* get_smatrix(string);

    void gnuplot_2D(unsigned int, unsigned int);
};

#endif
