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
#include <pthread.h>
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

#define uint unsigned int
#define HTK_FORMAT "htk"
#define XML_FORMAT "xml"
#define BORDER 10
#define COVMIN 1.0e-8
#define GAUSS_PUSH 0.1
#define TRANSMIN 1.0e-10

class ModelSet;
class Model;
class State;
class Stream;
class Gaussian;
class SVector;
class SMatrix;
class HMMLab_Object;
class StreamArea;
class FileData;

string gettag(istream&);
void init();
string execute(string cmd);
bool isallalpha(string);

struct point_len {
    uint i;
    uint j;
    double len;
    point_len(uint ii, uint jj, double l): i(ii), j(jj), len(l) {};
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

class FileData
{
    List<List<Vector*> > data;
    FileData(string, List<List<Vector*> > data);
public:
    bool selected;
    double maxprob;
    Model* model;
    string word;

    FileData();
    List<Vector*> operator[](uint);

    friend class Model;
    friend class ModelSet;
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
public:
    int index_distribution;
    double gconst;
    SVector* mean;
    SMatrix* covariance, *inv_covariance;
    List<uint> my_data;

    Gaussian(string, ModelSet*, int, double);
    Gaussian(string, ModelSet*, int, double, SVector*, SMatrix*);
    ~Gaussian();

    double probability(Vector*);
    void divide();
    void calc_gconst();

    friend class Stream;
    friend class ModelSet;
};


class Stream : public Shared
{
    List<Vector*> viterbi_data;
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
    double probability(Vector*);

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
    Gaussian* get_gaussian(uint, bool);
    bool has_gaussian(Gaussian*);

    bool in_viterbi_data(Vector*);
    void clear_viterbi_data();
    void add_viterbi_data(List<Vector*>);

    double probability(List<Vector*>);
    double probability_star(uint, List<Vector*>);

    friend class Model;
    friend class ModelSet;
};


class TransMatrix : public Shared
{
    void load(istream&, const char*, uint);
    void save(ostream&, const char*, uint);
    List<gsl_matrix*> matrix;
    List<TransMatrix*> joined_matrices;
public:
    TransMatrix(string, ModelSet*);
    TransMatrix(string, ModelSet*, int, double);
    ~TransMatrix();

    double operator()(uint, int);
    uint operator()(uint, int, double);
    void add_col_row();
    void remove_col_row(uint);
    TransMatrix* join_matrix(TransMatrix*);
    List<TransMatrix*> disjoint_matrix();

    friend class Model;
    friend class ModelSet;
};


class Model : public HMMLab_Object
{
    void load(istream&, const char*);
    void save(ostream&, const char*);
    List<Model*> joined_models;
public:
    ModelSet* modelset;
    List<State*> states;
    TransMatrix* trans_mat;

    Model(string, ModelSet*);
    Model(string, ModelSet*, List<State*>, TransMatrix*);
    ~Model();

    bool is_joined();
    Model* join_model(Model*);
    List<Model*> disjoint_model();

    void add_state(State*);
    void remove_state(State*);
    void remove_state(int);
    void select_gaussians();
    void unselect_gaussians();
    string create_image();
    void viterbi();
    void train(List<FileData*>);

    friend class ModelSet;
};

class StreamArea
{
    double screen_width, screen_height; //velkost DrawArea
    double graph_width, graph_height; //velkost grafu
    double graph_prob_width, graph_prob_height; //velkost grafu s pravdepodobnostami
    double pca_width, pca_height; //sirka a vyska PCA dat
    double edge_len_multiplier; //prisposoby dlzky hran medzi datami, aby najmensia dlzka bola 1

    List<double> edge_len; //vzdialenosti medzi datami

    List<Vector*> last_pos_data; //pozicie na grafe v poslednom layoute
    List<Vector*> last_gauss_pos; //pozicie stredov gaussianov v poslednom

    List<Vector*> last_pos_data_pca; //pozicie dat v poslednom vypocitanom PCA
    List<Vector*> last_gauss_pos_pca; //pozicie stredov gaussianov v poslednom vypocitanom PCA
    List<Vector*> last_gauss_var_pca; //variancie gaussianov v PCA

    List<Vector*> last_pos_data_prob; //pozicie dat v poslednom layoute s dlzkami hran podla pravdepodobnosti
    List<Vector*> last_gauss_pos_prob; //poslednom stredov gaussianov v poslednom layoute s dlzkami hran podla pravdepodobnosti

    List<Vector* >* get_positions(graph_t*, uint, const char*, bool prob);
    List<Vector* > translate_positions(List<Vector* >*);
    List<Vector* > translate_positions_prob(List<Vector* >*);
    List<Vector* > translate_pca_positions(List<Vector* >*, bool);
    graph_t* layout_graph(GVC_t*, bool);
    graph_t* layout_graph(GVC_t*, List<Vector*>);
    graph_t* layout_graph_prob(GVC_t*);
    Vector* get_pos(graph_t*, char*);

    void save_data_pos_2D(uint, string);
    void save_data_pos_3D(uint, uint, string);

public:
    ModelSet* modelset;

    List<Vector*> data; //data pripadajuce na tento stream

    List<Vector*> pos_data; //pozicie translatovane na velkost DrawArea
    List<Vector*> pos_gaussians;

    List<Vector*> pos_data_pca; //pozicie dat 2D PCA
    List<Vector*> pos_gaussians_pca; //pozicie stredov gaussianov 2D PCA
    List<Vector*> pos_gaussians_var_pca; //variancie gaussianov v 2D PCA

    List<Vector*> pos_data_prob;
    List<Vector*> pos_gaussians_prob;

    set<Gaussian*> selected_gaussians;

    StreamArea(ModelSet*);
    ~StreamArea();

    Vector* get_data(uint);

    void add_data(List<Vector*>);
    void set_wh(double, double);
    void reset_pos_gauss();
    void calc_pca();
    List<Vector*> get_data_2D(uint, uint);
    void calc_data_gauss();
    double calc_edge_len();

    friend class ModelSet;
};

class ModelSet : public HMMLab_Object
{
    void load(istream&, const char*);
    void save(ostream&, const char*);
    void create_cfg(string);
    void add_data(string, List<Vector*>);
public:
    uint dimension;
    uint streams_size;
    List<uint> streams_distribution;
    List<Model*> models;
    Dict<string, HMMLab_Object* > objects_dict;
    List<int> vecsize_tags;
    List<StreamArea*> stream_areas;
    set<Model*> drawarea_models;
    Dict<string, FileData*> files_data;

    ModelSet();
    ModelSet(string);
    ModelSet(string, const char*);
    ~ModelSet();
    void destroy();
    void save(const char*, const char*);

    void load_data(uint, string*);
    void reset_pos_gauss();

    void add_model(Model*);
    void remove_model(Model*);
    void remove_model(int);

    bool is_selected(Model*, int);
    bool is_selected(Gaussian*);
    uint selected_gaussians_count();
    uint loaded_data_count();
    List<Model*> get_models_with_gaussian(Gaussian*);
    string get_unique_name(string);

    bool gauss_cluster(List<Gaussian*>, List<Vector*>);
    bool gauss_push(bool, Gaussian*, Gaussian*);
    void train_model(Model*);

    void drawarea_models_append(Model*);
    void select_data(string);
    void unselect_data(string);

    Model* get_model(string);
    State* get_state(string);
    Stream* get_stream(string);
    Gaussian* get_gaussian(string);
    SVector* get_svector(string);
    SMatrix* get_smatrix(string);

    void gnuplot_2D(uint, uint);
};

#endif
