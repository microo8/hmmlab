#include "libhmm.h"


#ifndef HMMLABLIB_CPP
#define HMMLABLIB_CPP

/*-----------------Global-------------------*/
enum hmm_strings {
    default_s,
    global_options,
    streaminfo,
    vecsize,
    mfcc,
    mfcc_e,
    mfcc_e_d,
    mfcc_e_d_a,
    diagc,
    nulld,
    hmm_macro,
    begin_hmm,
    end_hmm,
    num_states,
    state,
    num_mixes,
    sweights,
    stream,
    mixture,
    mean_s,
    variance,
    gconst_s,
    transp,
    transp_macro,
    state_macro,
    variance_macro
};

static Dict<string, hmm_strings> hmm_strings_map;

void init()
{
    hmm_strings_map["default_s"] = default_s;
    hmm_strings_map["~o"] = global_options;
    hmm_strings_map["<STREAMINFO>"] = streaminfo;
    hmm_strings_map["<VECSIZE>"] = vecsize;
    hmm_strings_map["<MFCC>"] = mfcc;
    hmm_strings_map["<MFCC_E>"] = mfcc_e;
    hmm_strings_map["<MFCC_E_D>"] = mfcc_e_d;
    hmm_strings_map["<MFCC_E_D_A>"] = mfcc_e_d_a;
    hmm_strings_map["<DIAGC>"] = diagc;
    hmm_strings_map["<NULLD>"] = nulld;
    hmm_strings_map["~h"] = hmm_macro;
    hmm_strings_map["<BEGINHMM>"] = begin_hmm;
    hmm_strings_map["<ENDHMM>"] = end_hmm;
    hmm_strings_map["<NUMSTATES>"] = num_states;
    hmm_strings_map["<STATE>"] = state;
    hmm_strings_map["<NUMMIXES>"] = num_mixes;
    hmm_strings_map["<SWEIGHTS>"] = sweights;
    hmm_strings_map["<STREAM>"] = stream;
    hmm_strings_map["<MIXTURE>"] = mixture;
    hmm_strings_map["<MEAN>"] = mean_s;
    hmm_strings_map["<VARIANCE>"] = variance;
    hmm_strings_map["<GCONST>"] = gconst_s;
    hmm_strings_map["<TRANSP>"] = transp;
    hmm_strings_map["~t"] = transp_macro;
    hmm_strings_map["~s"] = state_macro;
    hmm_strings_map["~v"] = variance_macro;
};

class ValueException: public exception
{
    virtual const char* what() const throw() {
        return "Not good value";
    }
} ValEx;

string gettag(istream& in_stream)
{
    string result = "";
    char c;
    in_stream.get(c);
    while(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
        in_stream.get(c);
    }
    if(c != '<') {
        return "";
    }
    while(c != '>') {
        result.push_back(c);
        in_stream.get(c);
    }
    result.push_back(c);
    return result;
};

/*-----------------Global-------------------*/


/*-----------------Shared-------------------*/

Shared::Shared(string n, hmmlab_types t, ModelSet* ms): HMMLab_Object(n, t), modelset(ms) {};

void Shared::inc_ref_num()
{
    ref_num++;
    modelset->objects_dict[name] = this;
};

void Shared::dec_ref_num()
{
    ref_num--;
    if(ref_num <= 0) {
        modelset->objects_dict.erase(name);
        delete this;
    }
};

/*-----------------Shared-------------------*/


/*-----------------SVector------------------*/

SVector::SVector(string name, ModelSet* ms, int size, double elem = 0.0): Shared(name, SVECTOR, ms), Vector(size, elem) {};

void SVector::save(ostream& out_stream, const char* format)
{
    if(!strcmp(format, HTK_FORMAT)) {
        for(unsigned int i = 0; i < size(); i++) {
            out_stream << ' ' << scientific << (*this)[i];
        }
        out_stream << endl;
    } else if(!strcmp(format, XML_FORMAT)) {
    }
};
/*-----------------SVector------------------*/


/*-----------------SMatrix------------------*/

SMatrix::SMatrix(string name, ModelSet* ms, int m, int n, double elem = 0.0): Shared(name, SMATRIX, ms), Matrix(m, n, elem) {};

void SMatrix::save(ostream& out_stream, const char* format, bool all = false)
{
    if(!strcmp(format, HTK_FORMAT)) {
        if(all) {
            for(unsigned int i = 0; i < get_m(); i++) {
                for(unsigned int j = 0; j < get_n(); j++) {
                    out_stream << ' ' << scientific << (*this)(i, j);
                }
                out_stream << endl;
            }
        } else {
            for(unsigned int i = 0; i < get_m(); i++) {
                out_stream << ' ' << scientific << (*this)(i, i);
            }
            out_stream << endl;
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    }
};

/*-----------------SMatrix------------------*/


/*-----------------Gaussian-----------------*/

Gaussian::Gaussian(string name, ModelSet* ms, int index, double gc): Shared(name, GAUSSIAN, ms), index_distribution(index), gconst(gc) {};

Gaussian::Gaussian(string name, ModelSet* ms, int index, double gc, SVector* m, SMatrix* c): Shared(name, GAUSSIAN, ms), index_distribution(index), gconst(gc), mean(m), covariance(c)
{
    mean->inc_ref_num();
    covariance->inc_ref_num();
};

Gaussian::~Gaussian()
{
    mean->dec_ref_num();
    covariance->dec_ref_num();
};

void Gaussian::load(istream& in_stream, const char* format)
{
    bool gconst_readed = false;
    double x;
    int v;
    string line;
    char buffer[64];
    while(!in_stream.eof()) {
        in_stream >> skipws >> line;
        switch(hmm_strings_map[line]) {
        case mean_s:
            sprintf(buffer, "%s_mean", name.c_str());
            line = buffer;
            in_stream >> v;
            mean = new SVector(line, modelset, v, 0);
            for(int i = 0; i < v; i++) {
                in_stream >> scientific >> x;
                (*mean)(i, x);
            }
            break;
        case variance:
            sprintf(buffer, "%s_covariance", name.c_str());
            line = buffer;
            in_stream >> v;
            covariance = new SMatrix(line, modelset, v, v, 0);
            for(int i = 0; i < v; i++) {
                in_stream >> scientific >> x;
                (*covariance)(i, i, x);
            }
            break;
        case gconst_s:
            in_stream >> scientific >> gconst;
            gconst_readed = true;
            break;
        case variance_macro:
            in_stream >> skipws >> line;
            line = line.substr(1, line.length() - 2);
            covariance = static_cast<SMatrix*>(modelset->objects_dict[line]);
            covariance->inc_ref_num();
            break;
        default:
            in_stream.seekg((int)in_stream.tellg() - line.length(), ios::beg);
            line = "";
            break;
        }
        if(line == "") {
            break;
        }
    }
    if(!gconst_readed) {
        double sum = 0;
        int dim = modelset->streams_distribution[index_distribution];
        for(int i = 0; i < dim; i++) {
            sum += (*covariance)(i, i);
        }
        gconst = log(pow(2 * M_PI, dim) * sum);
    }
};

void Gaussian::save(ostream& out_stream, const char* format)
{
    if(!strcmp(format, HTK_FORMAT)) {
        out_stream << "<MEAN> " << mean->size() << endl;
        if(mean->ref_num > 1) {
            out_stream << "~u \"" << mean->name << '"' << endl;
        } else {
            mean->save(out_stream, HTK_FORMAT);
        }
        out_stream << "<VARIANCE> " << covariance->get_m() << endl;
        if(covariance->ref_num > 1) {
            out_stream << "~v \"" << covariance->name << '"' << endl;
        } else {
            covariance->save(out_stream, HTK_FORMAT);
        }
        out_stream << "<GCONST> " << gconst << endl;
    } else if(!strcmp(format, XML_FORMAT)) {
    }
};

/*-----------------Gaussian-----------------*/


/*------------------Stream------------------*/


Stream::Stream(string name, ModelSet* ms, int index): Shared(name, STREAM, ms), index_distribution(index), screen_width(0), screen_height(0) {};

Stream::Stream(string name, ModelSet* ms, int index, List<Gaussian*> g, List<double> g_w): Shared(name, STREAM, ms), index_distribution(index), screen_width(0), screen_height(0)
{
    if(g.size() != g_w.size()) {
        throw ValEx;
    }
    gaussians = g;
    gaussians_weights = g_w;
    for(unsigned int i = 0; i < gaussians.size(); i++) {
        gaussians[i]->inc_ref_num();
    }
};

Stream::~Stream()
{
    for(unsigned int i = 0; i < gaussians.size(); i++) {
        gaussians[i]->dec_ref_num();
    }
};

List<SVector*> Stream::gaussians_means()
{
    List<SVector*> result;
    for(unsigned int i = 0; i < gaussians.size(); i++) {
        result.append(gaussians[i]->mean);
    }
    return result;
};

void Stream::add_gaussian(Gaussian* g, double w)
{
    gaussians.append(g);
    gaussians_weights.append(w);
    g->inc_ref_num();
};

void Stream::remove_gaussian(Gaussian* g)
{
    int index = gaussians.index(g);
    remove_gaussian(index);
};

void Stream::remove_gaussian(int index)
{
    gaussians[index]->dec_ref_num();
    gaussians.remove(index);
    gaussians_weights.remove(index);
};

double Stream::get_gaussian_weight(Gaussian* g)
{
    return gaussians_weights[gaussians.index(g)];
};

void Stream::load(istream& in_stream, const char* format, int i_dist, int num_gaussians)
{
    unsigned int i;
    char buffer[64];
    double weight;
    string line;
    Gaussian* g;
    Dict<int, Gaussian*> gaussian_dict;
    Dict<int, double> gaussians_weights_dict;
    for(int j = 0; j < num_gaussians && !in_stream.eof(); j++) {
        in_stream >> skipws >> line;
        in_stream >> skipws >> i;
        in_stream >> skipws >> scientific >> weight;
        gaussians_weights_dict[i] = weight;
        sprintf(buffer, "%s_gaussian_%d", name.c_str(), i);
        line = buffer;
        g = new Gaussian(line, modelset, index_distribution, 0);
        g->load(in_stream, format);
        gaussian_dict[i] = g;
    }
    List<int> keys = gaussian_dict.keys();
    keys.listsort();
    for(i = 0; i < keys.size(); i++) {
        add_gaussian(gaussian_dict[keys[i]], gaussians_weights_dict[keys[i]]);
    }
};

graph_t* Stream::layout_graph(GVC_t* gvc)
{
    List<SVector* > gaussians_m = gaussians_means();

    char buffer [32];

    /* vytvori graf a prida kontrolne vrcholy */
    graph_t* g = agopen("", AGRAPHSTRICT);

    /* prida pozorovania a hrany medzi nimi */
    for(unsigned int i = 0; i < data.size(); i++) {
        sprintf(buffer, "node%d", i);
        Agnode_t* node = agnode(g, buffer);
        if(pos_data.size() > 0) {
            sprintf(buffer, "%8.6f,%8.6f", (*pos_data[i])[0], (*pos_data[i])[1]);
            agsafeset(node, "pin", "true", "");
            agsafeset(node, "pos", buffer, "");
        }
    }

    for(unsigned int i = 0; i < data.size(); i++) {
        int row = i / 2 * (2 * data.size() - i - 3) - 1;
        sprintf(buffer, "node%d", i);
        Agnode_t* i_node = agnode(g, buffer);
        for(unsigned int j = i + 1; j < data.size(); j++) {
            sprintf(buffer, "node%d", j);
            Agnode_t* j_node = agnode(g, buffer);
            Agedge_t* e = agedge(g, i_node, j_node);
            sprintf(buffer, "%8.6f", edge_len[row + j]);
            agsafeset(e, "len", buffer, "");
        }
    }

    /* prida stredy gaussianov a hrany medzi nimi a pozorovaniamy */
    for(unsigned int i = 0; i < gaussians_m.size(); i++) {
        sprintf(buffer, "gaussian%d", i);
        Agnode_t* gaussian = agnode(g, buffer);
        for(unsigned int j = 0; j < data.size(); j++) {
            sprintf(buffer, "node%d", j);
            Agnode_t* node = agnode(g, buffer);
            Agedge_t* e = agedge(g, gaussian, node);
            sprintf(buffer, "%8.6f", (*gaussians_m[i] - *data[j]).norm());
            agsafeset(e, "len", buffer, "");
        }
    }

    for(unsigned int i = 0; i < gaussians_m.size(); i++) {
        sprintf(buffer, "gaussian%d", i);
        Agnode_t* gaussian1 = agnode(g, buffer);
        for(unsigned int j = i + 1; j < gaussians_m.size(); j++) {
            sprintf(buffer, "gaussian%d", j);
            Agnode_t* gaussian2 = agnode(g, buffer);
            Agedge_t* e = agedge(g, gaussian1, gaussian2);
            sprintf(buffer, "%8.6f", (*gaussians_m[i] - *gaussians_m[j]).norm());
            agsafeset(e, "len", buffer, "");
        }
    }

    /* zavola neato na vypocitanie layoutu */
    gvLayout(gvc, g, "neato");
    attach_attrs(g);

    return g;
}

Vector* Stream::get_pos(graph_t* g, char* name_m)
{
    try {
        char* pos, *x, *y;
        int pch;

        Agnode_t* n = agnode(g, name_m);
        pos = agget(n, "pos");

        pch = strchr(pos, ',') - pos;
        x = new char[pch + 1];
        y = new char[strlen(pos) - pch + 1];

        strncpy(x, pos, pch);
        x[pch] = '\0';
        strncpy(y, pos + pch + 1, strlen(pos) - pch);
        y[strlen(pos) - pch] = '\0';

        Vector* result = new Vector(3, 1);
        (*result)(0, atof(x));
        (*result)(1, atof(y));
        delete x;
        delete y;
        return result;
    } catch(bad_alloc ex) {
        return NULL;
    }
}

List<Vector* > Stream::get_positions(graph_t* g, int size, const char* name)
{
    char buffer [16];

    char* bb =  agget(g, "bb") + 4;
    char* w, *h;
    int pch = strchr(bb, ',') - bb;
    w = new char[pch + 1];
    h = new char[strlen(bb) - pch + 1];
    strncpy(w, bb, pch);
    w[pch] = '\0';
    strncpy(h, bb + pch + 1, strlen(bb) - pch);
    h[strlen(bb) - pch] = '\0';
    double width = atof(w);
    double heigth = atof(h);
    delete[] w;
    delete[] h;

    double array [] = {screen_width / width, 0, 0, 0, screen_height / heigth, 0, 0, 0, 1};
    Matrix* translation = new Matrix(3, 3, 0);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            (*translation)(i, j, array[i * 3 + j]);
        }
    }
    List<Vector* > result(size, NULL);
    /* translacia kazdeho vrcholu */
    for(int i = 0; i < size; i++) {
        sprintf(buffer, "%s%d", name, i);
        Vector* vn = get_pos(g, buffer);
        if(vn != NULL) {
            result[i] = new Vector(3, 0);
            Vector tmp = (*translation) * (*vn);
            (*result[i])(0, tmp[0]);
            (*result[i])(1, tmp[1]);
            (*result[i])(2, tmp[2]);
        }
        delete vn;
    }
    delete translation;
    return result;
}

void Stream::add_data(List<Vector*> d)
{
    data = d;
    for(unsigned int i = 0; i < data.size(); i++) {
        unsigned int x = i / 2 * (2 * data.size() - i - 3);
        for(unsigned int j = i + 1; j < data.size(); j++) {
            edge_len[x + j - 1] = (*data[i] - *data[j]).norm();
        }
    }
    set_wh(screen_width, screen_height);
};

void Stream::set_wh(double w, double h)
{
    screen_width = w;
    screen_height = h;
    GVC_t* gvc = gvContext();
    graph_t* g = layout_graph(gvc);
    pos_data = get_positions(g, data.size(), "node");
    pos_gauss = get_positions(g, gaussians.size(), "gaussian");
    gvFreeLayout(gvc, g);
    agclose(g);
    gvFreeContext(gvc);
};

void Stream::save(ostream& out_stream, const char* format)
{
    if(!strcmp(format, HTK_FORMAT)) {
        if(gaussians.size() == 1) {
            //ak ma len jeden gaussian tak ho rovno ulozi
            if(gaussians[0]->ref_num > 1) {
                out_stream << "~m \"" << gaussians[0]->name << '"' << endl;
            } else {
                gaussians[0]->save(out_stream, HTK_FORMAT);
            }
        } else {
            for(unsigned int i = 0; i < gaussians.size(); i++) {
                out_stream << "<MIXTURE> " << i + 1 << ' ' << scientific << gaussians_weights[i] << endl;
                if(gaussians[i]->ref_num > 1) {
                    out_stream << "~m \"" << gaussians[i]->name << '"' << endl;
                } else {
                    gaussians[i]->save(out_stream, HTK_FORMAT);
                }
            }
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    }
};
/*------------------Stream------------------*/


/*-------------------State------------------*/

State::State(string name, ModelSet* ms): Shared(name, STATE, ms)
{
    for(int i = 0; i < ms->streams_size; i++) {
        char buffer[64];
        sprintf(buffer, "%s_stream_%d", name.c_str(), i);
        string name = buffer;
        Stream* s = new Stream(name, ms, i);
        s->inc_ref_num();
        streams.append(s);
        stream_weights.append(1.0 / ms->streams_size);
    }
};

State::State(string name, ModelSet* ms, List<double> s_w): Shared(name, STATE, ms)
{
    for(int i = 0; i < ms->streams_size; i++) {
        char buffer[64];
        sprintf(buffer, "%s_stream_%d", name.c_str(), i);
        string name = buffer;
        Stream* s = new Stream(name, ms, i);
        s->inc_ref_num();
        streams.append(s);
    }
    stream_weights = s_w;
};

State::State(string name, ModelSet* ms, List<Stream*> s, List<double> s_w): Shared(name, STATE, ms)
{
    if(s.size() != s_w.size()) {
        throw ValEx;
    }
    streams = s, stream_weights = s_w;
    for(unsigned int i = 0; i < streams.size(); i++) {
        streams[i]->inc_ref_num();
    }
};

State::~State()
{
    for(unsigned int i = 0; i < streams.size(); i++) {
        streams[i]->dec_ref_num();
    }
};

void State::load(istream& in_stream, const char* format)
{
    string line;
    int size;
    double weight;
    char buffer[128];
    Stream* str;
    List<int> num_gaussians_list;
    while(!in_stream.eof()) {
        in_stream >> skipws >> line;
        switch(hmm_strings_map[line]) {
        case num_mixes:
            for(int i = 0; i < modelset->streams_size; i++) {
                in_stream >> size;
                num_gaussians_list.append(size);
            }
            break;
        case sweights:
            in_stream >> skipws >> size;
            for(int i = 0; i < size; i++) {
                in_stream >> skipws >> weight;
                stream_weights[i] = weight;
            }
            break;
        case stream:
            in_stream >> skipws >> size;
            sprintf(buffer, "%s_stream_%d", name.c_str(), size);
            line = buffer;
            str = streams[size - 1];
            str->load(in_stream, format, size - 1, num_gaussians_list[size - 1]);
            break;
        default:
            if(hmm_strings_map[line] == mixture && modelset->streams_size == 1) {
                stream_weights.append(1.0);
                sprintf(buffer, "%s_stream_1", name.c_str());
                line = buffer;
                str = streams[0];
                str->load(in_stream, format, 0, num_gaussians_list[0]);
            } else {
                in_stream.seekg((int)in_stream.tellg() - line.length(), ios::beg);
            }
            line = "";
            break;
        }
        if(line == "") {
            break;
        }
    }
};

void State::save(ostream& out_stream, const char* format)
{
    if(!strcmp(format, HTK_FORMAT)) {
        out_stream << "<NUMMIXES>";
        unsigned int i;
        for(i = 0; i < streams.size(); i++) {
            out_stream << ' ' << streams[i]->gaussians.size();
        }
        out_stream << endl;
        if(streams.size() == 1) {
            //ak ma len jeden stream, tak ho ulozi bez
            //toho aby pre neho daval <STREAM> tag
            streams[0]->save(out_stream, HTK_FORMAT);
        } else {
            out_stream << "<SWEIGHTS> " << streams.size() << endl;
            for(i = 0; i < stream_weights.size(); i++) {
                out_stream << ' ' << scientific << stream_weights[i];
            }
            out_stream << endl;
            for(i = 0; i < streams.size(); i++) {
                out_stream << "<STREAM> " << i + 1 << endl;
                streams[i]->save(out_stream, HTK_FORMAT);
            }
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    }
};
/*-------------------State------------------*/


/*----------------TransMatrix---------------*/

TransMatrix::TransMatrix(string name, ModelSet* ms, int n, double value = 0): Shared(name, TRANSMATRIX, ms)
{
    List<List<double> * >* m = new List<List<double> * >();
    for(int i = 0; i < n; i++) {
        List<double>* l = new List<double>(n, value);
        m->append(l);
    }
    matrix.append(m);
};

TransMatrix::~TransMatrix()
{
    for(unsigned int i = 0; i < matrix.size(); i++) {
        for(unsigned int j = 0; j < matrix[i]->size(); j++) {
            delete(*matrix[i])[j];
        }
        delete matrix[i];
    }
};

double TransMatrix::operator()(unsigned int indexi, int indexj)
{
    unsigned int j = 0;
    for(j = 0; matrix[j]->size() <= indexi; j++) {
        if(j >= matrix.size()) {
            throw ValEx;
        }
        indexi -= matrix[j]->size();
        indexj -= matrix[j]->size();
    }
    if((unsigned int)indexj == matrix[j]->size() && indexi == matrix[j]->size() - 1) {
        return 1.0;
    }
    if(indexj < 0 || (unsigned int)indexj > matrix[j]->size()) {
        return 0.0;
    }
    return (*(*matrix[j])[indexi])[indexj];
};

void TransMatrix::operator()(unsigned int indexi, int indexj, double value)
{
    unsigned int j = 0;
    for(j = 0; matrix[j]->size() <= indexi; j++) {
        if(j >= matrix.size()) {
            throw ValEx;
        }
        indexi -= matrix[j]->size();
        indexj -= matrix[j]->size();
    }
    if((unsigned int)indexj == matrix[j]->size() && indexi == matrix[j]->size() - 1) {
        throw ValEx;
    }
    if(indexj < 0 || (unsigned int)indexj > matrix[j]->size()) {
        throw ValEx;
    }

    (*(*matrix[j])[indexi])[indexj] = value;
};

void TransMatrix::operator++()
{
    if(matrix.size() != 1) {
        throw ValEx;
    }
    for(unsigned int i = 0; i < matrix.size(); i++) {
        matrix[0][i].append(0);
    }
    matrix[0]->append(new List<double>(matrix.size() + 1, 0));
};

void TransMatrix::remove(int index)
{
    if(matrix.size() != 1) {
        throw ValEx;
    }
    matrix[0]->remove(index);
    for(unsigned int i = 0; i < matrix.size(); i++) {
        matrix[0][i].remove(index);
    }
};

void TransMatrix::remove_matrix(int index)
{
    matrix.remove(index);
};

void TransMatrix::add_matrix(TransMatrix& tm)
{
    for(unsigned int i = 0; i < tm.matrix.size(); i++) {
        matrix.append(tm.matrix[i]);
    }
};

void TransMatrix::load(istream& in_stream, const char* format, unsigned int size)
{
    double x;
    for(unsigned int j = 0; j < size; j++) {
        for(unsigned int k = 0; k < size; k++) {
            in_stream >> skipws >> scientific >> x;
            (*this)(j, k, x);
        }
    }
};

void TransMatrix::save(ostream& out_stream, const char* format, unsigned int size)
{
    if(!strcmp(format, HTK_FORMAT)) {
        for(unsigned int i = 0; i < size; i++) {
            for(unsigned int j = 0; j < size; j++) {
                out_stream << ' ' << scientific << (*this)(i, j);
            }
            out_stream << endl;
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    }
}
/*----------------TransMatrix---------------*/


/*-------------------Model------------------*/

Model::Model(string name, ModelSet* ms): HMMLab_Object(name, MODEL), modelset(ms)
{
    char buffer[64];
    sprintf(buffer, "%s_trans_matrix", name.c_str());
    string name_trans_mat = buffer;
    trans_mat = new TransMatrix(name_trans_mat, ms, 0);
    trans_mat->inc_ref_num();
};

Model::Model(string name, ModelSet* ms, List<State*> s, TransMatrix* t_m): HMMLab_Object(name, MODEL), modelset(ms), states(s), trans_mat(t_m)
{
    trans_mat->inc_ref_num();
    for(unsigned int i = 0; i < states.size(); i++) {
        states[i]->inc_ref_num();
    }
};

Model::~Model()
{
    for(unsigned int i = 0; i < states.size(); i++) {
        states[i]->dec_ref_num();
    }
    trans_mat->dec_ref_num();
};

void Model::add_state(State* s)
{
    states.append(s);
    s->inc_ref_num();
    trans_mat++;
};

void Model::remove_state(State* s)
{
    int index = states.index(s);
    remove_state(index);
};

void Model::remove_state(int index)
{
    states[index]->dec_ref_num();
    states.remove(index);
    trans_mat->remove(index);
};

void Model::load(istream& in_stream, const char* format)
{
    State* s;
    string line, state_name;
    int states_size = 0, i;
    char buffer[64];
    in_stream >> skipws >> line;
    name = line.substr(1, line.length() - 2);
    in_stream >> skipws >> line;
    Dict<int, State*> states_dict;
    while(hmm_strings_map[line] != end_hmm) {
        in_stream >> skipws >> line;
        switch(hmm_strings_map[line]) {
        case num_states:
            in_stream >> states_size;
            break;
        case state:
            in_stream >> i;
            in_stream >> skipws >> line;
            if(hmm_strings_map[line] == state_macro) {
                in_stream >> skipws >> line;
                line = line.substr(1, line.length() - 2);
                s = static_cast<State*>(modelset->objects_dict[line]);
                s->inc_ref_num();
            } else {
                in_stream.seekg((int)in_stream.tellg() - line.length(), ios::beg);
                sprintf(buffer, "%s_state_%d", name.c_str(), i);
                state_name = buffer;
                s = new State(state_name, modelset);
                s->load(in_stream, format);
            }
            states_dict[i] = s;
            break;
        case transp:
            sprintf(buffer, "%s_trans_mat", name.c_str());
            line = buffer;
            in_stream >> skipws >> i;
            trans_mat = new TransMatrix(line, modelset, i);
            trans_mat->load(in_stream, format, i);
            trans_mat->inc_ref_num();
            break;
        case transp_macro:
            in_stream >> skipws >> line;
            line = line.substr(1, line.length() - 2);
            trans_mat = static_cast<TransMatrix*>(modelset->objects_dict[line]);
            trans_mat->inc_ref_num();
            break;
        case end_hmm:
        default:
            line = "";
        }
        if(line == "") {
            break;
        }
    }
//vytvorit prvy a posledny state/////////////////////////////////////
    List<int> keys = states_dict.keys();
    keys.listsort();
    for(unsigned int i = 0; i < keys.size(); i++) {
        states.append(states_dict[keys[i]]);
    }
};

void Model::save(ostream& out_stream, const char* format)
{
    if(!strcmp(format, HTK_FORMAT)) {
        unsigned int states_size = states.size() + 2;
        out_stream << "<NUMSTATES> " << states_size << endl;
        for(unsigned int i = 0; i < states.size(); i++) {
            State* state = states[i];
            out_stream << "<STATE> " << i + 2 << endl;
            if(state->ref_num > 1) {
                out_stream << "~s \"" << state->name << '"' << endl;
            } else {
                state->save(out_stream, HTK_FORMAT);
            }
        }
        out_stream << "<TRANSP> " << states_size << endl;
        if(trans_mat->ref_num > 1) {
            out_stream << "~t \"" << trans_mat->name << '"' << endl;
        } else {
            trans_mat->save(out_stream, HTK_FORMAT, states_size);
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    }
};

/*-------------------Model------------------*/


/*------------------ModelSet----------------*/

ModelSet::ModelSet(): HMMLab_Object("modelset", MODELSET), streams_size(1) {};

ModelSet::ModelSet(string name): HMMLab_Object(name, MODELSET), streams_size(1) {};

ModelSet::ModelSet(string filename, const char* format): HMMLab_Object("modelset", MODELSET), streams_size(1)
{
    ifstream in_stream;
    in_stream.open(filename.c_str(), fstream::in);
    load(in_stream, format);
    in_stream.close();
};

ModelSet::~ModelSet()
{
    for(unsigned int i = 0; i < models.size(); i++) {
        delete models[i];
    }
};

void ModelSet::destroy()
{
    delete this;
};

void ModelSet::add_model(Model* m)
{
    objects_dict[m->name] = m;
    models.append(m);
};

void ModelSet::remove_model(Model* m)
{
    int index = models.index(m);
    remove_model(index);
};

void ModelSet::remove_model(int index)
{
    Model* m = models[index];
    objects_dict.erase(m->name);
    delete m;
    models.remove(index);
};

void ModelSet::load(istream& in_stream, const char* format)
{
    char c;
    int l;
    double x;
    string line, tmp_name;
    Model* model;
    TransMatrix* trans_mat;
    State* s;
    SMatrix* smat;

    while(!in_stream.eof()) {
        line = "";
        in_stream >> skipws >> line;
        switch(hmm_strings_map[line]) {
        case global_options:
            while(!in_stream.eof()) {
                in_stream >> skipws >> line;
                switch(hmm_strings_map[line]) {
                case streaminfo:
                    in_stream >> skipws >> streams_size;
                    dimension = 0;
                    for(int i = 0; i < streams_size; i++) {
                        int tmp_stream_size;
                        in_stream >> skipws >> tmp_stream_size;
                        dimension += tmp_stream_size;
                        streams_distribution.append(tmp_stream_size);
                    }
                    break;
                case vecsize:
                    in_stream.get(c);
                    while(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                        in_stream.get(c);
                    }
                    line = "";
                    while(c != ' ' && c != '\t' && c != '\n' && c != '\r' && c != '<') {
                        line.push_back(c);
                        in_stream.get(c);
                    }
                    if(streams_distribution.size() == 0) {
                        sscanf(line.c_str(), "%d", &l);
                        streams_distribution.append(l);
                    }
                    in_stream.seekg((int)in_stream.tellg() - 1, ios::beg);
                    while(!in_stream.eof()) {
                        line = gettag(in_stream);
                        switch(hmm_strings_map[line]) {
                        case mfcc:
                            vecsize_tags.append(mfcc);
                            break;
                        case mfcc_e_d_a:
                            vecsize_tags.append(mfcc_e_d_a);
                            break;
                        case nulld:
                            vecsize_tags.append(nulld);
                            break;
                        case diagc:
                            vecsize_tags.append(diagc);
                            break;
                        default:
                            in_stream.seekg((int)in_stream.tellg() - line.length() - 1, ios::beg);
                            line = "";
                            break;
                        }
                        if(line == "") {
                            break;
                        }
                    }
                    break;
                default:
                    in_stream.seekg((int)in_stream.tellg() - line.length(), ios::beg);
                    line == "";
                }
                if(line == "") {
                    break;
                }
            }
            line = "x";

            break;
        case hmm_macro:
            model = new Model("", this);
            model->load(in_stream, format);
            add_model(model);
            break;
        case transp_macro:
            in_stream >> skipws >> line;
            tmp_name = line.substr(1, line.length() - 2);
            in_stream >> skipws >> line;
            in_stream >> skipws >> l;
            trans_mat = new TransMatrix(tmp_name, this, l);
            trans_mat->load(in_stream, format, l);
            objects_dict[tmp_name] = trans_mat;
            break;
        case state_macro:
            in_stream >> skipws >> line;
            line = line.substr(1, line.length() - 2);
            s = new State(line, this);
            s->load(in_stream, format);
            objects_dict[line] = s;
            break;
        case variance_macro:
            in_stream >> skipws >> line;
            line = line.substr(1, line.length() - 2);
            in_stream >> skipws >> line;
            in_stream >> skipws >> l;
            smat = new SMatrix(line, this, l, l, 0);
            for(int i = 0; i < l; i++) {
                in_stream >> scientific >> x;
                (*smat)(i, i, x);
            }
            break;
        default:
            line == "";
        }
        if(line == "") {
            break;
        }
    }
};


void ModelSet::save(const char* filename, const char* format)
{
    ofstream out_stream;
    out_stream.open(filename, fstream::out);
    save(out_stream, format);
    out_stream.close();
};

void ModelSet::save(ostream& out_stream, const char* format)
{
    if(!strcmp(format, HTK_FORMAT)) {
        //najprv vlozi global data
        out_stream << "~o" << endl << "<STREAMINFO> " << streams_distribution.size() << ' ';
        for(unsigned int i = 0; i < streams_distribution.size(); i++) {
            out_stream << streams_distribution[i];
            if(i < streams_distribution.size() - 1) {
                out_stream << ' ';
            }
        }
        out_stream << endl << "<VECSIZE> " << sum(streams_distribution);
        for(unsigned int i = 0; i < vecsize_tags.size(); i++) {
            switch(vecsize_tags[i]) {
            case mfcc:
                out_stream << "<MFCC>";
                break;
            case mfcc_e:
                out_stream << "<MFCC_E>";
                break;
            case mfcc_e_d:
                out_stream << "<MFCC_E_D>";
                break;
            case mfcc_e_d_a:
                out_stream << "<MFCC_E_D_A>";
                break;
            case nulld:
                out_stream << "<NULLD>";
                break;
            case diagc:
                out_stream << "<DIAGC>";
                break;
            default:
                break;
            }
        }
        out_stream << endl;

        //prejde vsetky objekty a ak maju viac ako jeden pointer na seba
        //tak ich ulozi ako macro
        //zacina od najnizsich datovych struktur, keby niektore vyssie
        //potrebovali macro
        List<string> keys = objects_dict.keys();
        for(unsigned int i = 0; i < keys.size(); i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == SVECTOR && obj->ref_num > 1) {
                out_stream << "~u \"" << obj->name << '"' << endl;
                static_cast<SVector*>(obj)->save(out_stream, HTK_FORMAT);
            }
        }
        for(unsigned int i = 0; i < keys.size(); i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == SMATRIX && obj->ref_num > 1) {
                out_stream << "~v \"" << obj->name << '"' << endl;
                static_cast<SMatrix*>(obj)->save(out_stream, HTK_FORMAT);
            }
        }
        for(unsigned int i = 0; i < keys.size(); i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == TRANSMATRIX && obj->ref_num > 1) {
                out_stream << "~t \"" << obj->name << '"' << endl;
                unsigned int size = 0;
                for(unsigned int j = 0; j < keys.size(); j++) {
                    HMMLab_Object* om = objects_dict[keys[j]];
                    if(om->type == MODEL) {
                        Model* m = static_cast<Model*>(om);
                        if(m->trans_mat == obj) {
                            size = m->states.size();
                            break;
                        }
                    }
                }
                assert(size != 0);
                static_cast<TransMatrix*>(obj)->save(out_stream, HTK_FORMAT, size);
            }
        }
        for(unsigned int i = 0; i < keys.size(); i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == GAUSSIAN && obj->ref_num > 1) {
                out_stream << "~m \"" << obj->name << '"' << endl;
                static_cast<Gaussian*>(obj)->save(out_stream, HTK_FORMAT);
            }
        }
        for(unsigned int i = 0; i < keys.size(); i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == STATE && obj->ref_num > 1) {
                out_stream << "~s \"" << obj->name << '"' << endl;
                static_cast<State*>(obj)->save(out_stream, HTK_FORMAT);
            }
        }

        //uklada modely
        for(unsigned int i = 0; i < models.size(); i++) {
            Model* model = models[i];
            out_stream << "~h \"" << model->name << '"' << endl << "<BEGINHMM>" << endl;
            model->save(out_stream, HTK_FORMAT);
            out_stream << "<ENDHMM>" << endl;
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    };
};

Model* ModelSet::get_model(string name)
{
    return dynamic_cast<Model*>(objects_dict[name]);
};

State* ModelSet::get_state(string name)
{
    return dynamic_cast<State*>(objects_dict[name]);
};

Stream* ModelSet::get_stream(string name)
{
    return dynamic_cast<Stream*>(objects_dict[name]);
};

Gaussian* ModelSet::get_gaussian(string name)
{
    return dynamic_cast<Gaussian*>(objects_dict[name]);
};

SVector* ModelSet::get_svector(string name)
{
    return dynamic_cast<SVector*>(objects_dict[name]);
};

SMatrix* ModelSet::get_smatrix(string name)
{
    return dynamic_cast<SMatrix*>(objects_dict[name]);
};

/*------------------ModelSet----------------*/

#endif
