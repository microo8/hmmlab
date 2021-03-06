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

#include "libhmm.h"

#ifndef HMMLABLIB_CPP
#define HMMLABLIB_CPP

#define GRAPH_PROG "neato"

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

string execute(string cmd)
{
    FILE* pipe = popen(cmd.c_str(), "r");
    if(!pipe) {
        return "ERROR";
    }
    char buffer[128];
    string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL) {
            result += buffer;
        }
    }
    pclose(pipe);
    return result;
};

bool isallalpha(string str)
{
    for(uint i = 0; i < str.length(); i++) {
        if(!isalpha(str[i])) {
            return false;
        }
    }
    return true;
};

void* run_stream_area_pca(void* arg)
{
    ((StreamArea*)arg)->calc_pca();
    return NULL;
};

void* run_model_viterbi(void* arg)
{
    ((Model*)arg)->viterbi();
    return NULL;
}

template <typename T>
void run_threads(uint n, T* array, void * (*func)(void*))
{
    uint i = 0;
    void* status;
    pthread_t* thread = new pthread_t[n];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    for(i = 0; i < n; i++) {
        if(pthread_create(&thread[i], &attr, func, (void*)array[i])) {
            fprintf(stderr, "ERROR; return code from pthread_create()\n");
            exit(-1);
        }
    }
    pthread_attr_destroy(&attr);
    for(i = 0; i < n; i++) {
        if(pthread_join(thread[i], &status)) {
            fprintf(stderr, "ERROR; return code from pthread_join()\n");
            exit(-1);
        }
    }
    delete[] thread;
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
        for(uint i = 0; i < size(); i++) {
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
            for(uint i = 0; i < get_m(); i++) {
                for(uint j = 0; j < get_n(); j++) {
                    out_stream << ' ' << scientific << (*this)(i, j);
                }
                out_stream << endl;
            }
        } else {
            for(uint i = 0; i < get_m(); i++) {
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

void Gaussian::calc_gconst()
{
    double sum = 1;
    int dim = modelset->streams_distribution[index_distribution];
    for(int i = 0; i < dim; i++) {
        sum *= (*covariance)(i, i);
    }
    gconst = hmmlab_log(pow(2 * M_PI, dim) * sum);
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
            inv_covariance = new SMatrix(line + "_inv", modelset, v, v, 0);
            for(int i = 0; i < v; i++) {
                in_stream >> scientific >> x;
                (*covariance)(i, i, x);
                (*inv_covariance)(i, i, 1.0 / x);
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
            inv_covariance = static_cast<SMatrix*>(modelset->objects_dict[line + "_inv"]);
            inv_covariance->inc_ref_num();
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
        calc_gconst();
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

double Gaussian::probability(Vector* vec)
{
    assert(vec->size() == mean->size());
    double result;
    uint n = mean->size();
    gsl_vector* meang = mean->get_vector();
    gsl_vector* x = gsl_vector_alloc(n);
    gsl_vector_memcpy(x, vec->get_vector());
    gsl_matrix* winv = inv_covariance->get_matrix();
    gsl_vector_sub(x, meang);
    gsl_vector* tmp = gsl_vector_alloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1.0, winv, x, 0.0, tmp);
    gsl_blas_ddot(tmp, x, &result);
    gsl_vector_free(x);
    gsl_vector_free(tmp);
    result = (gconst + result) * -0.5;
    return result < EXPMIN ? EXPMIN : result;
};

void Gaussian::divide()
{
    string new_name;
    uint i, max_cov_index = -1;
    double max_cov = -DBL_MAX;
    List<Model*>::iterator mit;
    List<State*>::iterator sit;

    new_name = modelset->get_unique_name(name + "_copy");
    Gaussian* gauss = new Gaussian(new_name, modelset, index_distribution, gconst);
    modelset->objects_dict[new_name] = gauss;
    new_name = modelset->get_unique_name(mean->name + "_copy");
    gauss->mean = new SVector(new_name, modelset, mean->size(), 0.0);
    *gauss->mean = *mean;
    modelset->objects_dict[gauss->mean->name] = gauss->mean;
    new_name = modelset->get_unique_name(covariance->name + "_copy");
    gauss->covariance = new SMatrix(new_name, modelset, mean->size(), mean->size(), 0.0);
    *gauss->covariance = *covariance;
    modelset->objects_dict[gauss->covariance->name] = gauss->covariance;
    gauss->inv_covariance = new SMatrix(inv_covariance->name + "_copy", modelset, mean->size(), mean->size(), 0.0);
    *gauss->inv_covariance = *inv_covariance;
    for(i = 0; i < mean->size(); i++) {
        if((*covariance)(i, i) > max_cov) {
            max_cov = (*covariance)(i, i);
            max_cov_index = i;
        }
    }
    (*covariance)(max_cov_index, max_cov_index, max_cov + max_cov * 0.5);
    (*inv_covariance)(max_cov_index, max_cov_index, 1.0 / (max_cov + max_cov * 0.5));
    (*gauss->covariance)(max_cov_index, max_cov_index, max_cov - max_cov * 0.5);
    (*gauss->inv_covariance)(max_cov_index, max_cov_index, 1.0 / (max_cov - max_cov * 0.5));
    calc_gconst();
    gauss->calc_gconst();

    for(mit = modelset->models.begin(); mit < modelset->models.end(); mit++) {
        for(sit = (*mit)->states.begin(); sit < (*mit)->states.end(); sit++) {
            if((*sit)->streams[index_distribution]->gaussians.index(this) != -1) {
                (*sit)->streams[index_distribution]->add_gaussian(gauss, 1.0);
                modelset->stream_areas[index_distribution]->selected_gaussians.insert(gauss);
            }
        }
    }
    if(modelset->selected_gaussians_count() > 1 && modelset->loaded_data_count() > 0) {
        modelset->reset_pos_gauss();
    }
};
/*-----------------Gaussian-----------------*/


/*------------------Stream------------------*/


Stream::Stream(string name, ModelSet* ms, int index): Shared(name, STREAM, ms), index_distribution(index) {};

Stream::Stream(string name, ModelSet* ms, int index, List<Gaussian*> g, List<double> g_w): Shared(name, STREAM, ms), index_distribution(index)
{
    if(g.size() != g_w.size()) {
        throw ValEx;
    }
    gaussians = g;
    gaussians_weights = g_w;
    uint size = gaussians.size();
    for(uint i = 0; i < size; i++) {
        gaussians[i]->inc_ref_num();
    }
};

Stream::~Stream()
{
    uint size = gaussians.size();
    for(uint i = 0; i < size; i++) {
        gaussians[i]->dec_ref_num();
    }
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
    uint i;
    char buffer[64];
    double weight;
    string line;
    Gaussian* g;
    Dict<int, Gaussian*> gaussian_dict;
    Dict<int, double> gaussians_weights_dict;
    for(int j = 0; j < num_gaussians && !in_stream.eof(); j++) {
        in_stream >> skipws >> line;
        if(hmm_strings_map[line] == mixture) {
            in_stream >> skipws >> i;
            in_stream >> skipws >> scientific >> weight;
            gaussians_weights_dict[i] = weight;
            sprintf(buffer, "%s_gaussian_%d", name.c_str(), i);
            line = buffer;
            g = new Gaussian(line, modelset, index_distribution, 0);
            g->load(in_stream, format);
            gaussian_dict[i] = g;
        } else if(hmm_strings_map[line] == mean_s) {
            in_stream.seekg((int)in_stream.tellg() - line.length(), ios::beg);
            gaussians_weights_dict[i] = 1.0;
            sprintf(buffer, "%s_gaussian_%d", name.c_str(), 1);
            line = buffer;
            g = new Gaussian(line, modelset, index_distribution, 0);
            g->load(in_stream, format);
            gaussian_dict[i] = g;
        }
    }
    List<int> keys = gaussian_dict.keys();
    keys.listsort();
    uint size = keys.size();
    for(i = 0; i < size; i++) {
        add_gaussian(gaussian_dict[keys[i]], gaussians_weights_dict[keys[i]]);
    }
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
            uint size = gaussians.size();
            for(uint i = 0; i < size; i++) {
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

double Stream::probability(Vector* vec)
{
    uint i;
    double prob = 0.0, Gm = -DBL_MAX;
    uint size = gaussians.size();
    double* g = new double[size];
    for(i = 0; i < size; i++) {
        g[i] = hmmlab_log(gaussians_weights[i]) + gaussians[i]->probability(vec);
        if(Gm < g[i]) {
            Gm = g[i];
        }
    }
    for(i = 0; i < size; i++) {
        prob += hmmlab_exp(g[i] - Gm);
    }
    delete[] g;
    prob = Gm + hmmlab_log(prob);
    return prob < EXPMIN ? EXPMIN : prob;
};

/*------------------Stream------------------*/


/*-------------------State------------------*/

State::State(string name, ModelSet* ms): Shared(name, STATE, ms)
{
    for(uint i = 0; i < ms->streams_size; i++) {
        char buffer[64];
        sprintf(buffer, "%s_stream_%d", name.c_str(), i);
        string name = buffer;
        Stream* s = new Stream(name, ms, i);
        s->inc_ref_num();
        streams.append(s);
        stream_weights.append(1.0);
    }
};

State::State(string name, ModelSet* ms, List<double> s_w): Shared(name, STATE, ms)
{
    for(uint i = 0; i < ms->streams_size; i++) {
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
    uint size = streams.size();
    for(uint i = 0; i < size; i++) {
        streams[i]->inc_ref_num();
    }
};

State::~State()
{
    uint size = streams.size();
    for(uint i = 0; i < size; i++) {
        streams[i]->dec_ref_num();
    }
};

void State::select_gaussians(bool reset = false)
{
    uint size = streams.size();
    for(uint i = 0; i < size; i++) {
        StreamArea* strarea = modelset->stream_areas[i];
        List<Gaussian*>::iterator it;
        for(it = streams[i]->gaussians.begin(); it < streams[i]->gaussians.end(); it++) {
            strarea->selected_gaussians.insert(*it);
        };
        strarea->calc_data_gauss();
    }
    if(reset && modelset->selected_gaussians_count() > 1 && modelset->loaded_data_count() > 0) {
        modelset->reset_pos_gauss();
    }
};

void State::unselect_gaussians(bool reset = false)
{
    uint size = streams.size();
    for(uint i = 0; i < size; i++) {
        StreamArea* strarea = modelset->stream_areas[i];
        List<Gaussian*>::iterator it;
        for(it = streams[i]->gaussians.begin(); it < streams[i]->gaussians.end(); it++) {
            strarea->selected_gaussians.erase(*it);
        };
        strarea->calc_data_gauss();
    }
    if(reset && modelset->selected_gaussians_count() > 1 && modelset->loaded_data_count() > 1) {
        modelset->reset_pos_gauss();
    }
};

Gaussian* State::get_gaussian(uint index, bool select = false)
{
    List<Gaussian*>::iterator it;
    uint size = streams.size();
    for(uint i = 0; i < size; i++) {
        Stream* s = streams[i];
        for(it = s->gaussians.begin(); it < s->gaussians.end(); it++) {
            if(index-- == 0) {
                if(select) {
                    StreamArea* strarea = modelset->stream_areas[i];
                    if(strarea->selected_gaussians.count(*it) != 0) {
                        strarea->selected_gaussians.erase(*it);
                    } else {
                        strarea->selected_gaussians.insert(*it);
                    }
                    strarea->calc_data_gauss();
                    if(modelset->selected_gaussians_count() > 1 && modelset->loaded_data_count() > 1) {
                        modelset->reset_pos_gauss();
                    }
                }
                return *it;
            }
        }
    }
    return NULL;
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
            for(uint i = 0; i < modelset->streams_size; i++) {
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
                in_stream.seekg((int)in_stream.tellg() - line.length(), ios::beg);
                stream_weights.append(1.0);
                sprintf(buffer, "%s_stream_1", name.c_str());
                line = buffer;
                str = streams[0];
                str->load(in_stream, format, 0, num_gaussians_list[0]);
            } else if(hmm_strings_map[line] == mean_s) {
                in_stream.seekg((int)in_stream.tellg() - line.length(), ios::beg);
                stream_weights.append(1.0);
                sprintf(buffer, "%s_stream_1", name.c_str());
                line = buffer;
                str = streams[0];
                str->load(in_stream, format, 0, 1);
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
        uint i, size;
        size = streams.size();
        for(i = 0; i < size; i++) {
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

bool State::has_gaussian(Gaussian* g)
{
    List<Stream*>::iterator strit;
    for(strit = streams.begin(); strit < streams.end(); strit++) {
        if((*strit)->gaussians.index(g) != -1) {
            return true;
        }
    }
    return false;
};

bool State::in_viterbi_data(Vector* vec)
{
    List<Stream*>::iterator it;
    for(it = streams.begin(); it < streams.end(); it++) {
        if((*it)->viterbi_data.index(vec) != -1) {
            return true;
        }
    }
    return false;
};

void State::clear_viterbi_data()
{
    List<Stream*>::iterator it;
    for(it = streams.begin(); it < streams.end(); it++) {
        (*it)->viterbi_data.resize(0);
    }
};

void State::add_viterbi_data(List<Vector*> lvec)
{
    uint size = streams.size();
    for(uint i = 0; i < size; i++) {
        streams[i]->viterbi_data.append(lvec[i]);
    }
};

double State::probability(List<Vector*> lvec)
{
    assert(lvec.size() == streams.size());
    double prob = 0.0;
    uint size = streams.size();
    for(uint i = 0; i < size; i++) {
        prob += streams[i]->probability(lvec[i]) * stream_weights[i];
    }
    return prob < EXPMIN ? EXPMIN : prob;
};

double State::probability_star(uint stream_index, List<Vector*> lvec)
{
    assert(lvec.size() == streams.size());
    assert(stream_index < streams.size());
    double prob = 0.0;
    uint size = streams.size();
    for(uint i = 0; i < size; i++) {
        if(i != stream_index) {
            prob += streams[i]->probability(lvec[i]) * stream_weights[i];
        }
    }
    return prob < EXPMIN ? EXPMIN : prob;
};

/*-------------------State------------------*/


/*----------------TransMatrix---------------*/

TransMatrix::TransMatrix(string name, ModelSet* ms): Shared(name, TRANSMATRIX, ms) {};

TransMatrix::TransMatrix(string name, ModelSet* ms, int n, double value = 0): Shared(name, TRANSMATRIX, ms)
{
    gsl_matrix* m = gsl_matrix_alloc(n, n);
    gsl_matrix_set_all(m, value);
    matrix.append(m);
};

TransMatrix::~TransMatrix()
{
    if(joined_matrices.empty()) {
        for(uint i = 0; i < matrix.size(); i++) {
            gsl_matrix_free(matrix[i]);
        }
    }
};

double TransMatrix::operator()(uint indexi, int indexj)
{
    if(matrix.size() == 1) {
        return gsl_matrix_get(matrix[0], indexi, indexj);
    } else {
        uint j = 0;
        for(j = 0; (j < matrix.size()) && ((matrix[j]->size1 - (j == 0 ? 1 : 2)) <= indexi); j++) {
            indexi -= matrix[j]->size1 - (j == 0 ? 1 : 2);
            indexj -= matrix[j]->size1 - (j == 0 ? 1 : 2);
        }
        if(j >= matrix.size()) {
            return 0.0;
        }
        indexi += (j != 0);
        indexj += (j != 0);
        if(indexi < (uint)0 || indexj < 0 || (uint)matrix[j]->size1 <= indexi || matrix[j]->size1 <= (uint)indexj) {
            return 0.0;
        }
        return gsl_matrix_get(matrix[j], indexi, indexj);
    }
};

uint TransMatrix::operator()(uint indexi, int indexj, double value)
{
    if(matrix.size() == 1) {
        gsl_matrix_set(matrix[0], indexi, indexj, value);
        return 1;
    } else {
        uint j = 0;
        for(j = 0; (j < matrix.size()) && ((matrix[j]->size1 - (j == 0 ? 1 : 2)) <= indexi); j++) {
            indexi -= matrix[j]->size1 - (j == 0 ? 1 : 2);
            indexj -= matrix[j]->size1 - (j == 0 ? 1 : 2);
        }
        if(j >= matrix.size()) {
            return 0;
        }
        indexi += (j != 0);
        indexj += (j != 0);
        if(indexi < (uint)0 || indexj < 0 || (uint)matrix[j]->size1 <= indexi || matrix[j]->size1 <= (uint)indexj) {
            return 0;
        }
        gsl_matrix_set(matrix[j], indexi, indexj, value);
        return 1;
    }
};

void TransMatrix::add_col_row()
{
    assert(matrix.size() == 1);
    assert(ref_num == 1);
    uint size = matrix[0]->size1;
    gsl_matrix* m = gsl_matrix_calloc(size + 1, size + 1);
    gsl_matrix_view m1 = gsl_matrix_submatrix(m, 0, 0, size - 1, size - 1);
    gsl_matrix_view m2 = gsl_matrix_submatrix(matrix[0], 0, 0, size - 1, size - 1);
    gsl_matrix_memcpy(&m1.matrix, &m2.matrix);
    m1 = gsl_matrix_submatrix(m, 0, size, size - 1, 1);
    m2 = gsl_matrix_submatrix(matrix[0], 0, size - 1, size - 1, 1);
    gsl_matrix_memcpy(&m1.matrix, &m2.matrix);
    gsl_matrix_free(matrix[0]);
    matrix[0] = m;
};

void TransMatrix::remove_col_row(uint index)
{
    assert(matrix.size == 1);
    assert(ref_num == 1);
    assert(index < matrix[0]->size1);
    gsl_matrix* m1 = gsl_matrix_delete_row(matrix[0], index);
    gsl_matrix* m2 = gsl_matrix_delete_column(m1, index);
    gsl_matrix_free(matrix[0]);
    gsl_matrix_free(m1);
    matrix[0] = m2;
};

TransMatrix* TransMatrix::join_matrix(TransMatrix* tm)
{
    assert(tm->matrix.empty());
    modelset->objects_dict.erase(name);
    TransMatrix* ret = new TransMatrix(name + " " + tm->name, modelset);
    ret->matrix.resize(0);
    for(uint i = 0; i < matrix.size(); i++) {
        ret->matrix.append(matrix[i]);
    }
    ret->matrix.append(tm->matrix[0]);
    if(joined_matrices.empty()) {
        ret->joined_matrices.append(this);
    } else {
        ret->joined_matrices = joined_matrices;
        dec_ref_num();
    }
    modelset->objects_dict[ret->name] = ret;
    ret->joined_matrices.append(tm);
    return ret;
};

List<TransMatrix*> TransMatrix::disjoint_matrix()
{
    List<TransMatrix*> ret = joined_matrices;
    dec_ref_num();
    return ret;
};

void TransMatrix::load(istream& in_stream, const char* format, uint size)
{
    double x;
    for(uint j = 0; j < size; j++) {
        for(uint k = 0; k < size; k++) {
            in_stream >> skipws >> scientific >> x;
            (*this)(j, k, x);
        }
    }
};

void TransMatrix::save(ostream& out_stream, const char* format, uint size)
{
    if(!strcmp(format, HTK_FORMAT)) {
        for(uint i = 0; i < size; i++) {
            for(uint j = 0; j < size; j++) {
                out_stream << ' ' << scientific << (*this)(i, j);
            }
            out_stream << endl;
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    }
}
/*----------------TransMatrix---------------*/


/*-------------------Model------------------*/

Model::Model(string name, ModelSet* ms): HMMLab_Object(name, MODEL), modelset(ms), trans_mat(NULL) {};

Model::Model(string name, ModelSet* ms, List<State*> s, TransMatrix* t_m): HMMLab_Object(name, MODEL), modelset(ms), states(s), trans_mat(t_m)
{
    trans_mat->inc_ref_num();
    uint size = states.size();
    for(uint i = 0; i < size; i++) {
        states[i]->inc_ref_num();
    }
};

Model::~Model()
{
    uint size = states.size();
    for(uint i = 0; i < size; i++) {
        states[i]->dec_ref_num();
    }
    trans_mat->dec_ref_num();
};

void Model::add_state(State* s)
{
    states.append(s);
    s->inc_ref_num();
    trans_mat->add_col_row();
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
    trans_mat->remove_col_row(index + 1);
};

void Model::select_gaussians()
{
    List<State*>::iterator it;
    for(it = states.begin(); it < states.end(); it++) {
        (*it)->select_gaussians();
    }
};

void Model::unselect_gaussians()
{
    List<State*>::iterator it;
    for(it = states.begin(); it < states.end(); it++) {
        (*it)->unselect_gaussians();
    }
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
            s->inc_ref_num();
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
    List<int> keys = states_dict.keys();
    keys.listsort();
    uint size = keys.size();
    for(uint i = 0; i < size; i++) {
        states.append(states_dict[keys[i]]);
    }
};

void Model::save(ostream& out_stream, const char* format)
{
    if(!strcmp(format, HTK_FORMAT)) {
        uint states_size = states.size() + 2;
        out_stream << "<NUMSTATES> " << states_size << endl;
        uint size = states.size();
        for(uint i = 0; i < size; i++) {
            State* state = states[i];
            out_stream << "<STATE> " << i + 2 << endl;
            if(state->ref_num > 1) {
                out_stream << "~s \"" << state->name << '"' << endl;
            } else {
                state->save(out_stream, HTK_FORMAT);
            }
        }
        if(trans_mat->ref_num > 1) {
            out_stream << "~t \"" << trans_mat->name << '"' << endl;
        } else {
            out_stream << "<TRANSP> " << states_size << endl;
            trans_mat->save(out_stream, HTK_FORMAT, states_size);
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    }
};

string Model::create_image()
{
    GVC_t* gvc = gvContext();
    char n [256];
    memcpy(n, name.c_str(), name.size());

    graph_t* g = agopen(n, AGDIGRAPH);
    agsafeset(g, "rankdir", "LR", "");

    sprintf(n, "first_%s_state", name.c_str());
    Agnode_t* fnode = agnode(g, n);
    agsafeset(fnode, "style", "filled", "");
    agsafeset(fnode, "fillcolor", "gray", "");
    agsafeset(fnode, "filledsize", "true", "");
    agsafeset(fnode, "shape", "circle", "");
    memset(n, 0, 256);
    uint size = states.size();
    for(uint i = 0; i < size; i++) {
        sprintf(n, "%d", i + 1);
        Agnode_t* node = agnode(g, n);
        strcpy(n, states[i]->name.c_str());
        agsafeset(node, "label", n, "");
        agsafeset(node, "style", "filled", "");
        agsafeset(node, "fillcolor", "lightblue", "");
        agsafeset(node, "filledsize", "true", "");
        agsafeset(node, "shape", "circle", "");
    }
    memset(n, 0, 256);
    sprintf(n, "last_%s_state", name.c_str());
    Agnode_t* lnode = agnode(g, n);
    agsafeset(lnode, "style", "filled", "");
    agsafeset(lnode, "fillcolor", "gray", "");
    agsafeset(lnode, "filledsize", "true", "");
    agsafeset(lnode, "shape", "circle", "");

    memset(n, 0, 256);
    for(uint i = 0; i < size; i++) {
        if((*trans_mat)(0, i + 1) > 0) {
            sprintf(n, "%d", i + 1);
            Agnode_t* i_node = agnode(g, n);
            strcpy(n, states[i]->name.c_str());
            Agedge_t* e = agedge(g, fnode, i_node);
            sprintf(n, "%8.6f", (*trans_mat)(0, i + 1));
            agsafeset(e, "label", n, "");
        }
    }
    for(uint i = 0; i < size; i++) {
        if((*trans_mat)(i + 1, size + 1) > 0) {
            sprintf(n, "%d", i + 1);
            Agnode_t* i_node = agnode(g, n);
            Agedge_t* e = agedge(g, i_node, lnode);
            sprintf(n, "%8.6f", (*trans_mat)(i + 1, size + 1));
            agsafeset(e, "label", n, "");
        }
    }


    for(uint i = 0; i < size; i++) {
        for(uint j = 0; j < size; j++) {
            if((*trans_mat)(i + 1, j + 1) > 0) {
                sprintf(n, "%d", i + 1);
                Agnode_t* i_node = agnode(g, n);
                sprintf(n, "%d", j + 1);
                Agnode_t* j_node = agnode(g, n);
                Agedge_t* e = agedge(g, i_node, j_node);
                sprintf(n, "%8.6f", (*trans_mat)(i + 1, j + 1));
                agsafeset(e, "label", n, "");
            }
        }
    }

    gvLayout(gvc, g, "dot");
    //agwrite(g, stdout);
    char* path = tempnam("/dev/shm", NULL);
    FILE* fp = fopen(path, "wb");
    gvRender(gvc, g, "png", fp);
    fclose(fp);
    gvFreeLayout(gvc, g);
    agclose(g);
    gvFreeContext(gvc);
    string result = path;
    free(path);
    return result;
};

void Model::viterbi()
{
    uint i, j, k, t;
    List<State*>::iterator sit;
    map<string, FileData* >::iterator dit;
    List<Vector*> prob_list;

    for(sit = states.begin(); sit < states.end(); sit++) {
        (*sit)->clear_viterbi_data();
    }

    double maxprob, tmp;
    uint maxstate_index = 0;
    uint osize, size = states.size();
    double** b = new double*[size];
    uint** z = new uint*[size];
    double** psi = new double*[size];

    for(dit = modelset->files_data.begin(); dit != modelset->files_data.end(); dit++) {
        List<List<Vector*> > o = dit->second->data;
        osize = o[0].size();
        //vytvori b_i(o_j) -> log pravdepodobnost i-teho stavu a j-teho vektora
        for(i = 0; i < size; i++) {
            b[i] = new double[osize];
            z[i] = new uint[osize];
            psi[i] = new double[osize];
            for(j = 0; j < osize; j++) {
                prob_list.resize(0);
                for(k = 0; k < o.size(); k++) {
                    prob_list.append(o[k][j]);
                }
                b[i][j] = states[i]->probability(prob_list);
            }
        }

        //prvy krok \psi_1(1) = 1
        psi[0][0] = 0;
        //druhy prok \psi_j(1) = a_{1j}b_j(o_1) -> log(a_{1j}) + log(b_j(o_1))
        for(j = 1; j < size; j++) {
            psi[j][0] = hmmlab_log((*trans_mat)(0, j)) + b[j][0];
        }
        //treti krok \psi_j(t) = \max_i{\psi_i(t-1) + log(a_{ij})} + log(b_j(o_t))
        for(t = 1; t < osize; t++) {
            for(j = 0; j < size; j++) {
                maxprob = -DBL_MAX;
                for(i = 0; i < size; i++) {
                    tmp = psi[i][t - 1] + hmmlab_log((*trans_mat)(i + 1, j + 1));
                    if(maxprob < tmp) {
                        maxprob = tmp;
                        maxstate_index = i;
                    }
                }
                psi[j][t] = maxprob + b[j][t];
                z[j][t] = maxstate_index;
            }
        }
        maxprob = -DBL_MAX;
        for(i = 0; i < size; i++) {
            tmp = psi[i][osize - 1];
            if(maxprob < tmp) {
                maxprob = tmp;
                maxstate_index = i;
            }
        }
        prob_list.resize(0);
        for(k = 0; k < o.size(); k++) {
            prob_list.append(o[k][osize - 1]);
        }
        states[maxstate_index]->add_viterbi_data(prob_list);
        for(i = osize - 2; i > 0; i--) {
            prob_list.resize(0);
            for(k = 0; k < o.size(); k++) {
                prob_list.append(o[k][i]);
            }
            states[maxstate_index = z[maxstate_index][i]]->add_viterbi_data(prob_list);
        }

        for(i = 0; i < size; i++) {
            delete[] b[i];
            delete[] z[i];
            delete[] psi[i];
        }
        if(dit->second->model == NULL || dit->second->maxprob < maxprob) {
            dit->second->model = this;
            dit->second->maxprob = maxprob;
        }
    }

    delete[] b;
    delete[] z;
    delete[] psi;
};

bool Model::is_joined()
{
    return !joined_models.empty();
};

Model* Model::join_model(Model* m)
{
    assert(m.joined_models.empty());
    if(joined_models.empty()) {
        List<State*> joined_states = states;
        joined_states += m->states;
        Model* joined_model = new Model(name + "_" + m->name, modelset, joined_states, trans_mat->join_matrix(m->trans_mat));
        modelset->objects_dict[joined_model->name] = joined_model;
        joined_model->inc_ref_num();
        joined_model->joined_models.append(this);
        joined_model->joined_models.append(m);
        modelset->drawarea_models_append(joined_model);
        modelset->drawarea_models.erase(this);
        modelset->drawarea_models.erase(m);
        modelset->reset_pos_gauss();
        return joined_model;
    } else {
        states += m->states;
        modelset->objects_dict.erase(name);
        name = name + "_" + m->name;
        modelset->objects_dict[name] = this;
        trans_mat = trans_mat->join_matrix(m->trans_mat);
        joined_models.append(m);
        modelset->drawarea_models.erase(m);
        modelset->reset_pos_gauss();
        return this;
    }
};

List<Model*> Model::disjoint_model()
{
    List<Model*> ret = joined_models;
    modelset->drawarea_models.erase(this);
    dec_ref_num();
    return ret;
};

void Model::train(List<FileData*> data)
{
    double prob;
    int t;
    uint i, j, r;
    uint N = states.size() + 2;
    uint R = data.size();
    gsl_vector* d = gsl_vector_alloc(N - 2);
    gsl_vector* tmpvec = gsl_vector_alloc(R);

    gsl_matrix* a = gsl_matrix_alloc(N, N);
    for(i = 0; i < N; i++) {
        for(j = 0; j < N; j++) {
            gsl_matrix_set(a, i, j, hmmlab_log((*trans_mat)(i, j)));
        }
    }

    gsl_matrix** alpha = new gsl_matrix*[R];
    gsl_matrix** beta = new gsl_matrix*[R];
    gsl_matrix** U = new gsl_matrix*[R];
    double* P = new double[R];
    uint* T = new uint[R];

    for(r = 0; r < R; r++) {
        T[r] = data[r]->data[0].size();

        alpha[r] = gsl_matrix_calloc(N, T[r]);

        //prvy krok \alpha_1(1) = 1
        gsl_matrix_set(alpha[r], 0, 0, 0);

        //druhy krok \alpha_j(1) = a_{1j}*b_j(o_1)
        for(j = 1; j < N - 1; j++) {
            gsl_matrix_set(alpha[r], j, 0, gsl_matrix_get(a, 0, j) + states[j - 1]->probability((*data[r])[0]));
        }

        //pre 1 < j < N a 1 < t < T treti krok:
        //\alpha_j(t) = [\sum_{i=2}^{N-1} \aplha_i(t-1)*a_{ij}] * b_j(o_t)
        for(t = 1; t < T[r]; t++) {
            for(j = 1; j < N - 1; j++) {
                for(i = 1; i < N - 1; i++) {
                    gsl_vector_set(d, i - 1, gsl_matrix_get(alpha[r], i, t - 1) + gsl_matrix_get(a, i, j));
                }
                gsl_matrix_set(alpha[r], j, t, logsumexp(d) + states[j - 1]->probability((*data[r])[t]));
            }
        }

        //posledny krok \alpha_N(T) = \sum_{i=2}^{N-1} \alpha_i(T) * a_{iN}
        for(i = 1; i < N - 1; i++) {
            gsl_vector_set(d, i - 1, gsl_matrix_get(alpha[r], i, T[r] - 1) + gsl_matrix_get(a, i, N - 1));
        }
        gsl_matrix_set(alpha[r], N - 1, T[r] - 1, logsumexp(d));

        beta[r] = gsl_matrix_calloc(N, T[r]);
        //prvy krok \beta_i(T) = a_{iN}
        for(i = 0; i < N; i++) {
            gsl_matrix_set(beta[r], i, T[r] - 1, gsl_matrix_get(a, i, N - 1));
        }
        //druhy krok pre 1 < i < N a T > t >= 1 je \beta_i(t) = \sum_{j=2}^{N-1} a_{ij} * b_j(o_{t+1})\beta_j(t+1)
        for(t = T[r] - 2; t >= 0; --t) {
            for(i = 1; i < N - 1; i++) {
                for(j = 1; j < N - 1; j++) {
                    gsl_vector_set(d, j - 1, gsl_matrix_get(a, i, j) + states[j - 1]->probability((*data[r])[t + 1]) + gsl_matrix_get(beta[r], j, t + 1));
                }
                gsl_matrix_set(beta[r], i, t, logsumexp(d));
            }
        }
        //poledny krok \beta_1(1) = \sum_{j=2}^{N-1} a_{1j} * b_j(o_1) * \beta_j(1)
        for(j = 1; j < N - 1; j++) {
            gsl_vector_set(d, j - 1, gsl_matrix_get(a, 0, j) + states[j - 1]->probability((*data[r])[0]) + gsl_matrix_get(beta[r], j, 0));
        }
        P[r] = logsumexp(d);
        gsl_matrix_set(beta[r], 0, 0, P[r]);
    }

    //prestavenie prechodovych pravdepodobnosti
    for(i = 0; i < trans_mat->matrix.size(); i++) {
        gsl_matrix_set_all(trans_mat->matrix[i], 0.0);
    }
    //prvy krok a_{1j} = \frac{1}{R} \sum_{r=1}^{R} \frac{1}{P_r} \alpha_j(1) * \beta_j(1)
    for(j = 1; j < N - 1; j++) {
        for(r = 0; r < R; r++) {
            gsl_vector_set(tmpvec, r, -P[r] + gsl_matrix_get(alpha[r], j, 0) + gsl_matrix_get(beta[r], j, 0));
        }
        prob = hmmlab_exp(logsumexp(tmpvec) - hmmlab_log(R));
        (*trans_mat)(0, j, prob < TRANSMIN ? 0.0 : prob);
    }

    //druhy krok
    for(i = 1; i < N - 1; i++) {
        for(r = 0; r < R; r++) {
            gsl_vector* tmpsum = gsl_vector_alloc(T[r]);
            for(t = 0; t < T[r]; t++) {
                gsl_vector_set(tmpsum, t, gsl_matrix_get(alpha[r], i, t) + gsl_matrix_get(beta[r], i, t));
            }
            gsl_vector_set(tmpvec, r, -P[r] + logsumexp(tmpsum));
            gsl_vector_free(tmpsum);
        }
        double under = logsumexp(tmpvec);
        for(j = 1; j < N - 1; j++) {
            for(r = 0; r < R; r++) {
                gsl_vector* tmpsum = gsl_vector_alloc(T[r] - 1);
                for(t = 0; t < T[r] - 1; t++) {
                    gsl_vector_set(tmpsum, t, gsl_matrix_get(alpha[r], i, t) +
                                   gsl_matrix_get(a, i, j) +
                                   states[j - 1]->probability((*data[r])[t + 1]) +
                                   gsl_matrix_get(beta[r], j, t + 1));
                }
                gsl_vector_set(tmpvec, r, -P[r] + logsumexp(tmpsum));
                gsl_vector_free(tmpsum);
            }
            prob = hmmlab_exp(logsumexp(tmpvec) - under);
            (*trans_mat)(i, j, prob < TRANSMIN ? 0.0 : prob);
        }
    }

    //posledny
    for(i = 1; i < N - 1; i++) {
        for(r = 0; r < R; r++) {
            gsl_vector_set(tmpvec, r, -P[r] + gsl_matrix_get(alpha[r], i, T[r] - 1) + gsl_matrix_get(beta[r], i, T[r] - 1));
        }
        prob = logsumexp(tmpvec);
        for(r = 0; r < R; r++) {
            gsl_vector* tmpsum = gsl_vector_alloc(T[r]);
            for(t = 0; t < T[r]; t++) {
                gsl_vector_set(tmpsum, t, gsl_matrix_get(alpha[r], i, t) + gsl_matrix_get(beta[r], i, t));
            }
            gsl_vector_set(tmpvec, r, -P[r] + logsumexp(tmpsum));
            gsl_vector_free(tmpsum);
        }
        prob = hmmlab_exp(prob - logsumexp(tmpvec));
        (*trans_mat)(i, N - 1, prob < TRANSMIN ? 0.0 : prob);
    }

    for(r = 0; r < R; r++) {
        U[r] = gsl_matrix_alloc(N - 2, T[r]);
        for(j = 0; j < N - 2; j++) {
            for(t = 0; t < T[r]; t++) {
                if(t == 0) {
                    gsl_matrix_set(U[r], j, 0, gsl_matrix_get(a, 0, j + 1));
                } else {
                    for(i = 1; i < N - 1; i++) {
                        gsl_vector_set(d, i - 1, gsl_matrix_get(alpha[r], i, t - 1) + gsl_matrix_get(a, i, j + 1));
                    }
                    gsl_matrix_set(U[r], j, t, logsumexp(d));
                }
            }
        }
    }

    double**** *L = new double**** [r];
    for(r = 0; r < R; r++) {
        L[r] = new double** *[N - 2];
        for(j = 0; j < N - 2; j++) {
            L[r][j] = new double** [modelset->streams_size];
            for(uint s = 0; s < modelset->streams_size; s++) {
                L[r][j][s] = new double*[states[j]->streams[s]->gaussians.size()];
                for(uint m = 0; m < states[j]->streams[s]->gaussians.size(); m++) {
                    L[r][j][s][m] = new double[T[r]];
                    double c_jsm = hmmlab_log(states[j]->streams[s]->gaussians_weights[m]);
                    for(t = 0; t < T[r]; t++) {
                        L[r][j][s][m][t] = -P[r] +
                                           gsl_matrix_get(U[r], j, t) +
                                           c_jsm +
                                           states[j]->streams[s]->gaussians[m]->probability((*data[r])[t][s]) +
                                           gsl_matrix_get(beta[r], j + 1, t) +
                                           states[j]->probability_star(s, (*data[r])[t]);
                    }
                }
            }
        }
    }

    //prepocita nove stredy vsetkych gaussianov
    for(j = 0; j < N - 2; j++) {
        for(uint s = 0; s < modelset->streams_size; s++) {
            gsl_vector* tmp_data = gsl_vector_alloc(modelset->streams_distribution[s]);
            for(uint m = 0; m < states[j]->streams[s]->gaussians.size(); m++) {
                gsl_vector* mean = states[j]->streams[s]->gaussians[m]->mean->get_vector();
                gsl_vector_set_all(mean, 0.0);
                double Lsum = 0.0;
                for(r = 0; r < R; r++) {
                    for(t = 0; t < T[r]; t++) {
                        gsl_vector_memcpy(tmp_data, (*data[r])[t][s]->get_vector());
                        gsl_vector_scale(tmp_data, hmmlab_exp(L[r][j][s][m][t]));
                        gsl_vector_add(mean, tmp_data);
                        Lsum += hmmlab_exp(L[r][j][s][m][t]);
                    }
                }
                gsl_vector_scale(mean, 1.0 / Lsum);
            }
            gsl_vector_free(tmp_data);
        }
    }

    //prepocita kovariancne matice vsetkych gaussianov
    for(j = 0; j < N - 2; j++) {
        for(uint s = 0; s < modelset->streams_size; s++) {
            gsl_vector* tmp_data = gsl_vector_alloc(modelset->streams_distribution[s]);
            for(uint m = 0; m < states[j]->streams[s]->gaussians.size(); m++) {
                double Lsum = 0.0;
                gsl_matrix* cov = states[j]->streams[s]->gaussians[m]->covariance->get_matrix();
                gsl_matrix_set_all(cov, 0.0);
                for(r = 0; r < R; r++) {
                    for(t = 0; t < T[r]; t++) {
                        gsl_vector_memcpy(tmp_data, (*data[r])[t][s]->get_vector());
                        gsl_vector_sub(tmp_data, states[j]->streams[s]->gaussians[m]->mean->get_vector());
                        gsl_blas_dger(hmmlab_exp(L[r][j][s][m][t]), tmp_data, tmp_data, cov);
                        Lsum += hmmlab_exp(L[r][j][s][m][t]);
                    }
                }
                gsl_matrix_scale(cov, 1.0 / Lsum);
                //upravi na diagonalnu
                for(uint ii = 0; ii < cov->size1; ii++) {
                    for(uint jj = 0; jj < cov->size2; jj++) {
                        if(ii == jj) {
                            prob = gsl_matrix_get(cov, ii, ii);
                            gsl_matrix_set(cov, ii, ii, prob < COVMIN ? COVMIN : prob);
                        } else {
                            gsl_matrix_set(cov, ii, jj, 0.0);
                        }
                    }
                }
                //prepocita inverznu
                gsl_matrix* tmpmat = states[j]->streams[s]->gaussians[m]->inv_covariance->get_matrix();
                for(uint jj = 0; jj < cov->size1; jj++) {
                    gsl_matrix_set(tmpmat, jj, jj, 1.0 / gsl_matrix_get(cov, jj, jj));
                }
                //prepocita GCONST
                states[j]->streams[s]->gaussians[m]->calc_gconst();
            }
            gsl_vector_free(tmp_data);
        }
    }

    //prepocita vahy gaussianov
    for(j = 0; j < N - 2; j++) {
        double under = 0.0;
        for(uint s = 0; s < modelset->streams_size; s++) {
            for(uint m = 0; m < states[j]->streams[s]->gaussians.size(); m++) {
                for(r = 0; r < R; r++) {
                    for(t = 0; t < T[r]; t++) {
                        under += hmmlab_exp(L[r][j][s][m][t]);
                    }
                }
            }
        }
        for(uint s = 0; s < modelset->streams_size; s++) {
            for(uint m = 0; m < states[j]->streams[s]->gaussians.size(); m++) {
                double c_jsm = 0.0;
                for(r = 0; r < R; r++) {
                    for(t = 0; t < T[r]; t++) {
                        c_jsm += hmmlab_exp(L[r][j][s][m][t]);
                    }
                }
                states[j]->streams[s]->gaussians_weights[m] = c_jsm / under;
            }
        }
    }

    for(r = 0; r < R; r++) {
        for(j = 0; j < N - 2; j++) {
            for(uint s = 0; s < modelset->streams_size; s++) {
                for(uint m = 0; m < states[j]->streams[s]->gaussians.size(); m++) {
                    delete[] L[r][j][s][m];
                }
                delete[] L[r][j][s];
            }
            delete[] L[r][j];
        }
        delete[] L[r];
    }

    gsl_matrix_free(a);
    gsl_vector_free(d);
    gsl_vector_free(tmpvec);
    for(r = 0; r < R; r++) {
        gsl_matrix_free(alpha[r]);
        gsl_matrix_free(beta[r]);
        gsl_matrix_free(U[r]);
    }
    delete[] alpha;
    delete[] beta;
    delete[] U;
    delete[] P;
    delete[] T;
    delete[] L;
};

/*-------------------Model------------------*/

/*----------------StreamArea----------------*/

StreamArea::StreamArea(ModelSet* ms, bool g)
{
    modelset = ms;
    selected_gaussians.clear();
    graphviz = g;
};

StreamArea::~StreamArea()
{
    /*
    data
    last_pos_data
    last_pos_data_pca
    last_gauss_pos
    last_gauss_pos_pca
    last_gauss_var_pca
    last_pos_data_prob
    last_gauss_pos_prob
    pos_data
    pos_data_pca
    pos_gaussians
    pos_gaussians_pca
    pos_gaussians_var_pca
    pos_data_prob
    pos_gaussians_prob
    */
    uint size = data.size();
    for(uint i = 0; i < size; i++) {
        delete data[i];
    }
    data.resize(0);

    size = last_pos_data.size();
    for(uint i = 0; i < size; i++) {
        delete last_pos_data[i];
    }
    last_pos_data.resize(0);
    size = last_pos_data_pca.size();
    for(uint i = 0; i < size; i++) {
        delete last_pos_data_pca[i];
    }
    last_pos_data_pca.resize(0);

    size = last_gauss_pos.size();
    for(uint i = 0; i < size; i++) {
        delete last_gauss_pos[i];
    }
    last_gauss_pos.resize(0);

    size = last_gauss_pos_pca.size();
    for(uint i = 0; i < size; i++) {
        delete last_gauss_pos_pca[i];
    }
    last_gauss_pos_pca.resize(0);

    size = last_gauss_var_pca.size();
    for(uint i = 0; i < size; i++) {
        delete last_gauss_var_pca[i];
    }
    last_gauss_var_pca.resize(0);

    size = last_pos_data_prob.size();
    for(uint i = 0; i < size; i++) {
        delete last_pos_data_prob[i];
    }
    last_pos_data_prob.resize(0);

    size = last_gauss_pos_prob.size();
    for(uint i = 0; i < size; i++) {
        delete last_gauss_pos_prob[i];
    }
    last_gauss_pos_prob.resize(0);

    size = pos_data.size();
    for(uint i = 0; i < size; i++) {
        delete pos_data[i];
    }
    pos_data.resize(0);

    size = pos_data_pca.size();
    for(uint i = 0; i < size; i++) {
        delete pos_data_pca[i];
    }
    pos_data_pca.resize(0);

    size = pos_gaussians.size();
    for(uint i = 0; i < size; i++) {
        delete pos_gaussians[i];
    }
    pos_gaussians.resize(0);

    size = pos_gaussians_pca.size();
    for(uint i = 0; i < size; i++) {
        delete pos_gaussians_pca[i];
    }
    pos_gaussians_pca.resize(0);

    size = pos_gaussians_var_pca.size();
    for(uint i = 0; i < size; i++) {
        delete pos_gaussians_var_pca[i];
    }
    pos_gaussians_var_pca.resize(0);

    size = pos_data_prob.size();
    for(uint i = 0; i < size; i++) {
        delete pos_data_prob[i];
    }
    pos_data_prob.resize(0);

    size = pos_gaussians_prob.size();
    for(uint i = 0; i < size; i++) {
        delete pos_gaussians_prob[i];
    }
    pos_gaussians_prob.resize(0);
};

void StreamArea::set_wh(double w, double h)
{
    screen_width = w;
    screen_height = h;
    uint size;

    //vycisti pos_data
    size = pos_data.size();
    for(uint i = 0; i < size; i++) {
        delete pos_data[i];
    }
    pos_data.resize(0);
    pos_data = translate_positions(&last_pos_data);

    size = pos_gaussians.size();
    for(uint i = 0; i < size; i++) {
        delete pos_gaussians[i];
    }
    pos_gaussians.resize(0);
    pos_gaussians = translate_positions(&last_gauss_pos);

    size = pos_data_pca.size();
    for(uint i = 0; i < size; i++) {
        delete pos_data_pca[i];
    }
    pos_data_pca.resize(0);
    pos_data_pca = translate_pca_positions(&last_pos_data_pca, true);

    size = pos_gaussians_pca.size();
    for(uint i = 0; i < size; i++) {
        delete pos_gaussians_pca[i];
    }
    pos_gaussians_pca.resize(0);
    pos_gaussians_pca = translate_pca_positions(&last_gauss_pos_pca, true);

    size = pos_gaussians_var_pca.size();
    for(uint i = 0; i < size; i++) {
        delete pos_gaussians_var_pca[i];
    }
    pos_gaussians_var_pca.resize(0);
    pos_gaussians_var_pca = translate_pca_positions(&last_gauss_var_pca, false);

    size = pos_data_prob.size();
    for(uint i = 0; i < size; i++) {
        delete pos_data_prob[i];
    }
    pos_data_prob.resize(0);
    pos_data_prob = translate_positions_prob(&last_pos_data_prob);

    size = pos_gaussians_prob.size();
    for(uint i = 0; i < size; i++) {
        delete pos_gaussians_prob[i];
    }
    pos_gaussians_prob.resize(0);
    pos_gaussians_prob = translate_positions_prob(&last_gauss_pos_prob);
};

graph_t* StreamArea::layout_graph(GVC_t* gvc, bool run = false)
{
    if(graphviz) {
        char buffer [256];
        double elen;
        uint size;

        /* vytvori graf a prida kontrolne vrcholy */
        graph_t* g = agopen("", AGRAPHSTRICT);
        agsafeset(g, "overlap", "scale", "");

        size = data.size();
        for(uint i = 0; i < size; i++) {
            int row = i / 2 * (2 * size - i - 3) - 1;
            sprintf(buffer, "node%d", i);
            Agnode_t* i_node = agnode(g, buffer);
            for(uint j = i + 1; j < size; j++) {
                sprintf(buffer, "node%d", j);
                Agnode_t* j_node = agnode(g, buffer);
                Agedge_t* e = agedge(g, i_node, j_node);
                elen = edge_len[row + j] * edge_len_multiplier;
                sprintf(buffer, "%8.6f", elen == 0 ? 1 : elen);
                agsafeset(e, "len", buffer, "");
            }
        }

        if(run) {
            gvLayout(gvc, g, GRAPH_PROG);
            attach_attrs(g);
            //agwrite(g, stdout);
        }

        return g;
    } else {
        return NULL;
    }
};

graph_t* StreamArea::layout_graph(GVC_t* gvc, List<Vector* > gaussians_m)
{
    if(graphviz) {
        uint size, size2;
        char buffer [256];
        graph_t* g = layout_graph(gvc);
        double minimum = 1, elen;
        List<struct point_len> gauss_data;
        size = gaussians_m.size();
        for(uint i = 0; i < size; i++) {
            size2 = data.size();
            for(uint j = 0; j < size2; j++) {
                double len = (*gaussians_m[i] - *data[j]).norm();
                struct point_len t(i, j, len);
                gauss_data.append(t);
                if(len < 1 && len > 0  && minimum > len) {
                    minimum = len;
                }
            }
        }

        List<struct point_len> gauss_gauss;
        size = gaussians_m.size();
        for(uint i = 0; i < size; i++) {
            for(uint j = i + 1; j < size; j++) {
                double len = (*gaussians_m[i] - *gaussians_m[j]).norm();
                struct point_len t(i, j, len);
                gauss_gauss.append(t);
                if(len < 1 && len > 0  && minimum > len) {
                    minimum = len;
                }
            }
        }

        if(edge_len_multiplier < (1 / minimum)) {
            edge_len_multiplier = 1 / minimum;
        }
        /* prida stredy gaussianov a hrany medzi nimi a pozorovaniamy */
        size = gauss_data.size();
        for(uint index = 0; index < size; index++) {
            uint i = gauss_data[index].i;
            uint j = gauss_data[index].j;
            sprintf(buffer, "gaussian%d", i);
            Agnode_t* gaussian = agnode(g, buffer);
            sprintf(buffer, "node%d", j);
            Agnode_t* node = agnode(g, buffer);
            Agedge_t* e = agedge(g, gaussian, node);
            elen = gauss_data[index].len * edge_len_multiplier;
            sprintf(buffer, "%8.6f", elen == 0 ? 1 : elen);
            agsafeset(e, "len", buffer, "");
        }

        size = gauss_gauss.size();
        for(uint index = 0; index < size; index++) {
            uint i = gauss_gauss[index].i;
            uint j = gauss_gauss[index].j;
            sprintf(buffer, "gaussian%d", i);
            Agnode_t* gaussian1 = agnode(g, buffer);
            sprintf(buffer, "gaussian%d", j);
            Agnode_t* gaussian2 = agnode(g, buffer);
            Agedge_t* e = agedge(g, gaussian1, gaussian2);
            elen = gauss_gauss[index].len * edge_len_multiplier;
            sprintf(buffer, "%8.6f", elen == 0 ? 1 : elen);
            agsafeset(e, "len", buffer, "");
        }

        /* zavola neato na vypocitanie layoutu */
        gvLayout(gvc, g, GRAPH_PROG);
        attach_attrs(g);
        //agwrite(g, stdout);

        return g;
    } else {
        return NULL;
    }
}

graph_t* StreamArea::layout_graph_prob(GVC_t* gvc)
{
    if(graphviz) {
        uint size, size2;
        char buffer [256];
        graph_t* g = agopen("", AGRAPHSTRICT);
        agsafeset(g, "overlap", "scale", "");
        List<Gaussian*> list_selected_gaussians = set2List(selected_gaussians);
        if(list_selected_gaussians.size() == 0) {
            gvFreeLayout(gvc, g);
            agclose(g);
            return NULL;
        }
        double elen, prob, minprob = DBL_MAX, maxprob = -DBL_MAX;
        List<struct point_len> gauss_data;
        size = list_selected_gaussians.size();
        for(uint i = 0; i < size; i++) {
            Gaussian* g = list_selected_gaussians[i];
            size2 = data.size();
            for(uint j = 0; j < size2; j++) {
                prob = -g->probability(data[j]);
                if(minprob > prob) {
                    minprob = prob;
                }
                if(maxprob < prob) {
                    maxprob = prob;
                }
                struct point_len t(i, j, prob);
                gauss_data.append(t);
            }
        }
        List<struct point_len> gauss_gauss;
        size = list_selected_gaussians.size();
        for(uint i = 0; i < size; i++) {
            Gaussian* gauss1 = list_selected_gaussians[i];
            for(uint j = i + 1; j < size; j++) {
                Gaussian* gauss2 = list_selected_gaussians[j];
                double prob1 = gauss1->probability(gauss2->mean);
                double prob2 = gauss2->probability(gauss1->mean);
                double prob = -(prob1 > prob2 ? prob1 : prob2);
                if(minprob > prob) {
                    minprob = prob;
                }
                if(maxprob < prob) {
                    maxprob = prob;
                }
                struct point_len t(i, j, prob);
                gauss_gauss.append(t);
            }
        }

        size = gauss_data.size();
        for(uint i = 0; i < size; i++) {
            gauss_data[i].len -= minprob;
            gauss_data[i].len *= 100000.0 / (maxprob - minprob);
            gauss_data[i].len += 1.0;
        }
        size = gauss_gauss.size();
        for(uint i = 0; i < size; i++) {
            gauss_gauss[i].len -= minprob;
            gauss_gauss[i].len *= 100000.0 / (maxprob - minprob);
            gauss_gauss[i].len += 1.0;
        }
        /* prida stredy gaussianov a hrany medzi nimi a pozorovaniamy */
        size = gauss_data.size();
        for(uint index = 0; index < size; index++) {
            uint i = gauss_data[index].i;
            uint j = gauss_data[index].j;
            sprintf(buffer, "gaussian%d", i);
            Agnode_t* gaussian = agnode(g, buffer);
            sprintf(buffer, "node%d", j);
            Agnode_t* node = agnode(g, buffer);
            Agedge_t* e = agedge(g, gaussian, node);
            elen = gauss_data[index].len;
            sprintf(buffer, "%8.6f", elen == 0 ? 1 : elen);
            agsafeset(e, "len", buffer, "");
        }

        size = gauss_gauss.size();
        for(uint index = 0; index < size; index++) {
            uint i = gauss_gauss[index].i;
            uint j = gauss_gauss[index].j;
            sprintf(buffer, "gaussian%d", i);
            Agnode_t* gaussian1 = agnode(g, buffer);
            sprintf(buffer, "gaussian%d", j);
            Agnode_t* gaussian2 = agnode(g, buffer);
            Agedge_t* e = agedge(g, gaussian1, gaussian2);
            elen = gauss_gauss[index].len;
            sprintf(buffer, "%8.6f", elen == 0 ? 1 : elen);
            agsafeset(e, "len", buffer, "");
        }

        /* zavola neato na vypocitanie layoutu */
        gvLayout(gvc, g, GRAPH_PROG);
        attach_attrs(g);
        //agwrite(g, stdout);

        return g;
    } else {
        return NULL;
    }
}

Vector* StreamArea::get_pos(graph_t* g, char* name_m)
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

        Vector* result = new Vector(3, 0);

        stringstream ssx(stringstream::in | stringstream::out);
        stringstream ssy(stringstream::in | stringstream::out);
        ssx << x;
        ssy << y;
        double pos_x;
        ssx >> pos_x;
        double pos_y;
        ssy >> pos_y;
        (*result)(0, pos_x);
        (*result)(1, pos_y);
        delete x;
        delete y;
        return result;
    } catch(bad_alloc ex) {
        return NULL;
    }
}

List<Vector* >* StreamArea::get_positions(graph_t* g, uint size, const char* name, bool prob)
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
    stringstream ssw(stringstream::in | stringstream::out);
    ssw << scientific << w;
    if(prob) {
        ssw >> scientific >> graph_prob_width;
    } else {
        ssw >> scientific >> graph_width;
    }
    stringstream ssh(stringstream::in | stringstream::out);
    ssh << scientific << h;
    if(prob) {
        ssh >> scientific >> graph_prob_height;
    } else {
        ssh >> scientific >> graph_height;
    }
    delete[] w;
    delete[] h;

    List<Vector* >* result = new List<Vector*>(size, NULL);
    for(uint i = 0; i < size; i++) {
        sprintf(buffer, "%s%d", name, i);
        (*result)[i] = get_pos(g, buffer);
    }
    return result;
};

List<Vector*> StreamArea::translate_positions(List<Vector*>* veclist)
{
    uint size;
    double array [] = {(screen_width - BORDER * 2) / graph_width, 0, 0, 0, (screen_height - BORDER * 2) / graph_height, 0, 0, 0, 1};
    Matrix translation(3, 3, 0);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            translation(i, j, array[i * 3 + j]);
        }
    }
    size = veclist->size();
    List<Vector* > result(size, NULL);
    /* translacia kazdeho vrcholu */
    for(uint i = 0; i < size; i++) {
        result[i] = new Vector(3, 0);
        Vector tmp = translation * *(*veclist)[i];
        (*result[i])(0, tmp[0] + BORDER);
        (*result[i])(1, tmp[1] + BORDER);
        (*result[i])(2, tmp[2]);
    }
    return result;
};

List<Vector*> StreamArea::translate_positions_prob(List<Vector*>* veclist)
{
    uint size;
    double array [] = {(screen_width - BORDER * 2) / graph_prob_width, 0, 0, 0, (screen_height - BORDER * 2) / graph_prob_height, 0, 0, 0, 1};
    Matrix translation(3, 3, 0);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            translation(i, j, array[i * 3 + j]);
        }
    }
    size = veclist->size();
    List<Vector* > result(size, NULL);
    /* translacia kazdeho vrcholu */
    for(uint i = 0; i < size; i++) {
        result[i] = new Vector(3, 0);
        Vector tmp = translation * *(*veclist)[i];
        (*result[i])(0, tmp[0] + BORDER);
        (*result[i])(1, tmp[1] + BORDER);
        (*result[i])(2, tmp[2]);
    }
    return result;
};


List<Vector*> StreamArea::translate_pca_positions(List<Vector*>* veclist, bool center)
{
    uint size;
    double array [] = {(screen_width - BORDER * 2) / pca_width, 0, 0, 0, (screen_height - BORDER * 2) / pca_height, 0, 0, 0, 1};
    Matrix translation(3, 3, 0);
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            translation(i, j, array[i * 3 + j]);
        }
    }
    Vector screen_center(3, 0);
    screen_center(0, screen_width / 2.0);
    screen_center(1, screen_height / 2.0);
    size = veclist->size();
    List<Vector* > result(size, NULL);
    /* translacia kazdeho vrcholu */
    for(uint i = 0; i < size; i++) {
        result[i] = new Vector(3, 0);
        Vector tmp = translation * *(*veclist)[i];
        if(center) {
            tmp += screen_center;
        }
        (*result[i])(0, tmp[0] + BORDER);
        (*result[i])(1, tmp[1] + BORDER);
        (*result[i])(2, tmp[2]);
    }
    return result;
};

double StreamArea::calc_edge_len()
{
    uint size;
    double edge_len_minimum = 1.0;
    //vymaze vypocitane dlzky hran dat a vypocita nove
    edge_len.resize(0);
    size = data.size();
    for(uint i = 0; i < size; i++) {
        for(uint j = i + 1; j < size; j++) {
            Vector sub = *data[i] - *data[j];
            double len = sub.norm();
            edge_len.append(len);
            if(len < 1 && len > 0 && edge_len_minimum > len) {
                edge_len_minimum = len;
            }
        }
    }
    edge_len_multiplier = 1.0 / edge_len_minimum;
    return edge_len_minimum;
};

void StreamArea::add_data(List<Vector*> d)
{
    uint size;
    data += d;
    calc_edge_len();

    //vymaze posledne pozicie
    size = last_pos_data.size();
    for(uint i = 0; i < size; i++) {
        delete last_pos_data[i];
    }
    last_pos_data.resize(0);
    size = last_pos_data_prob.size();
    for(uint i = 0; i < size; i++) {
        delete last_pos_data_prob[i];
    }
    last_pos_data_prob.resize(0);
    if(graphviz) {
        if(selected_gaussians.size() == 0) {
            //graphvizom vypocita pozicie dlzoveho grafu
            void* status;
            pthread_t thread;
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            if(pthread_create(&thread, &attr, run_stream_area_pca, (void*)this)) {
                fprintf(stderr, "ERROR; return code from pthread_create()\n");
                exit(-1);
            }
            pthread_attr_destroy(&attr);

            GVC_t* gvc = gvContext();
            graph_t* g = layout_graph(gvc, true);
            List<Vector*>* list = get_positions(g, data.size(), "node", false);
            last_pos_data = *list;
            delete list;
            gvFreeLayout(gvc, g);
            agclose(g);
            gvFreeContext(gvc);

            //graphvizom vypocita pozicie pravdepodobnostneho grafu
            size = last_pos_data_prob.size();
            for(uint i = 0; i < size; i++) {
                delete last_pos_data_prob[i];
            }
            last_pos_data_prob.resize(0);
            set_wh(screen_width, screen_height);

            if(pthread_join(thread, &status)) {
                fprintf(stderr, "ERROR; return code from pthread_join()\n");
                exit(-1);
            }
        } else {
            reset_pos_gauss();
        }
    }
};


void StreamArea::reset_pos_gauss()
{
    void* status;
    pthread_t thread;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    if(pthread_create(&thread, &attr, run_stream_area_pca, (void*)this)) {
        fprintf(stderr, "ERROR; return code from pthread_create()\n");
        exit(-1);
    }
    pthread_attr_destroy(&attr);

    if(graphviz) {
        GVC_t* gvc = gvContext();
        set<Gaussian*>::iterator itg;
        List<Vector*>::iterator it;
        List<Vector*> gauss_means;
        for(itg = selected_gaussians.begin(); itg != selected_gaussians.end(); itg++) {
            gauss_means.append((*itg)->mean);
        }
        graph_t* g = layout_graph(gvc, gauss_means);
        //gaussian pos
        List<Vector*>* list = get_positions(g, gauss_means.size(), "gaussian", false);
        for(it = last_gauss_pos.begin(); it < last_gauss_pos.end(); it++) {
            delete *it;
            *it = NULL;
        }
        last_gauss_pos.resize(0);
        last_gauss_pos = *list;

        //gauss pos del and trans
        for(it = pos_gaussians.begin(); it < pos_gaussians.end(); it++) {
            delete *it;
            *it = NULL;
        }
        pos_gaussians.resize(0);
        pos_gaussians = translate_positions(list);
        delete list;

        //node pos
        list = get_positions(g, data.size(), "node", false);
        for(it = last_pos_data.begin(); it < last_pos_data.end(); it++) {
            delete *it;
            *it = NULL;
        }
        last_pos_data.resize(0);
        last_pos_data = *list;

        //node pos del and trans
        for(it = pos_data.begin(); it < pos_data.end(); it++) {
            delete *it;
            *it = NULL;
        }
        pos_data.resize(0);
        pos_data = translate_positions(list);
        delete list;

        gvFreeLayout(gvc, g);
        agclose(g);
        gvFreeContext(gvc);

        gvc = gvContext();
        g = layout_graph_prob(gvc);
        //gaussian pos prob
        for(it = last_gauss_pos_prob.begin(); it < last_gauss_pos_prob.end(); it++) {
            delete *it;
            *it = NULL;
        }
        last_gauss_pos_prob.resize(0);
        //gauss pos prob del and trans
        for(it = pos_gaussians_prob.begin(); it < pos_gaussians_prob.end(); it++) {
            delete *it;
            *it = NULL;
        }
        pos_gaussians_prob.resize(0);
        if(g != NULL) {
            list = get_positions(g, gauss_means.size(), "gaussian", true);
            last_gauss_pos_prob = *list;
            pos_gaussians_prob = translate_positions_prob(list);
            delete list;
        }

        //node pos
        for(it = last_pos_data_prob.begin(); it < last_pos_data_prob.end(); it++) {
            delete *it;
            *it = NULL;
        }
        last_pos_data_prob.resize(0);
        //node pos del and trans
        for(it = pos_data_prob.begin(); it < pos_data_prob.end(); it++) {
            delete *it;
            *it = NULL;
        }
        pos_data_prob.resize(0);
        if(g != NULL) {
            list = get_positions(g, data.size(), "node", true);
            last_pos_data_prob = *list;
            pos_data_prob = translate_positions_prob(list);
            delete list;

            gvFreeLayout(gvc, g);
            agclose(g);
        }
        gvFreeContext(gvc);
        if(pthread_join(thread, &status)) {
            fprintf(stderr, "ERROR; return code from pthread_join()\n");
            exit(-1);
        }
    }
};

void StreamArea::save_data_pos_2D(uint dim, string filename)
{
    assert(dim < modelset->streams_distribution[modelset->stream_areas.index(this)]);
    uint size;
    double x, y, sigma, mu;
    double flipedpi = 1.0 / sqrt(2 * M_PI);
    set<Gaussian*>::iterator it;
    ofstream data_file(("/dev/shm/" + filename).c_str());

    size = data.size();
    for(uint i = 0; i < size; i++) {
        x = (*data[i])[dim];
        data_file << scientific << x << ' ';

        y = 0;
        for(it = selected_gaussians.begin(); it != selected_gaussians.end(); it++) {
            mu = (*((*it)->mean))[dim];
            sigma = (*((*it)->covariance))(dim, dim);
            //1./(sigma*sqrt(2*pi)) * exp( -(x-mu)**2 / (2*sigma**2) )
            y  += flipedpi * (1.0 / sqrt(sigma)) * hmmlab_exp(-pow(x - mu, 2) / (2 * sigma));
        }
        data_file << scientific << y << endl;
        data_file << scientific << x << ' ' << scientific << 0 << endl;
    }

    data_file.close();
};

List<Vector*> StreamArea::get_data_2D(uint dim1, uint dim2)
{
    List<Vector*> result;
    List<Vector*>::iterator it;
    for(it = data.begin(); it < data.end(); it++) {
        Vector* v = *it;
        Vector* n = new Vector(2);
        (*n)(0, (*v)[dim1]);
        (*n)(1, (*v)[dim2]);
        result.append(n);
    }
    return result;
};

void StreamArea::calc_pca()
{
    uint M = modelset->streams_distribution[modelset->stream_areas.index(this)];
    uint N = data.size();
    List<Vector*>::iterator it;
    set<Gaussian*>::iterator git;
    uint i = 0;
    if(N > 0) {

        //vytvori maticu
        gsl_matrix* m = gsl_matrix_alloc(M, N);
        for(it = data.begin(); it < data.end(); it++) {
            Vector* v = *it;
            gsl_vector* vv = v->get_vector();
            gsl_matrix_set_col(m, i++, vv);
        }

        //spusti pca
        gsl_matrix* pca_m = gsl_pca(m, 2);
        //vycisti doterajsie data
        for(it = last_pos_data_pca.begin(); it < last_pos_data_pca.end(); it++) {
            delete *it;
            *it = NULL;
        }
        last_pos_data_pca.resize(0);

        for(it = last_gauss_pos_pca.begin(); it < last_gauss_pos_pca.end(); it++) {
            delete *it;
            *it = NULL;
        }
        last_gauss_pos_pca.resize(0);

        for(it = last_gauss_var_pca.begin(); it < last_gauss_var_pca.end(); it++) {
            delete *it;
            *it = NULL;
        }
        last_gauss_var_pca.resize(0);

        //vynasobi novu bazu s maticou dat
        gsl_matrix* data_pca = gsl_matrix_alloc(2, N);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, pca_m, m, 0.0, data_pca);
        gsl_matrix_free(m);

        //vypocita vysku, sirku a priemer
        gsl_vector_view rx = gsl_matrix_row(data_pca, 0);
        pca_width = gsl_vector_max(&rx.vector) - gsl_vector_min(&rx.vector);
        gsl_vector_view ry = gsl_matrix_row(data_pca, 1);
        pca_height = gsl_vector_max(&ry.vector) - gsl_vector_min(&ry.vector);
        gsl_vector* pca_mean = gsl_vector_alloc(2);
        gsl_vector_set(pca_mean, 0, (gsl_vector_max(&rx.vector) + gsl_vector_min(&rx.vector)) / 2.0);
        gsl_vector_set(pca_mean, 1, (gsl_vector_max(&ry.vector) + gsl_vector_min(&ry.vector)) / 2.0);
        for(i = 0; i < N; i++) {
            gsl_vector_view data_col = gsl_matrix_column(data_pca, i);
            gsl_vector_sub(&data_col.vector, pca_mean);
        }

        //priradi nove data
        uint size = data.size();
        for(i = 0; i < size; i++) {
            Vector* v = new Vector(3, 0);
            (*v)(0, gsl_matrix_get(data_pca, 0, i));
            (*v)(1, gsl_matrix_get(data_pca, 1, i));
            last_pos_data_pca.append(v);
        }
        for(git = selected_gaussians.begin(); git != selected_gaussians.end(); git++) {
            gsl_vector* gv = (*git)->mean->get_vector();
            gsl_vector* ggv = gsl_vector_alloc(2);
            gsl_blas_dgemv(CblasTrans, 1.0, pca_m, gv, 0.0, ggv);
            gsl_vector_sub(ggv, pca_mean);
            Vector* v = new Vector(3, 0);
            (*v)(0, gsl_vector_get(ggv, 0));
            (*v)(1, gsl_vector_get(ggv, 1));
            last_gauss_pos_pca.append(v);
            gsl_vector_free(ggv);
        }
        gsl_vector_free(pca_mean);
        gsl_matrix_free(data_pca);

        //prepocita variancie
        for(git = selected_gaussians.begin(); git != selected_gaussians.end(); git++) {
            gsl_vector* diag_covariance = gsl_vector_alloc(M);
            Gaussian* gauss = *git;
            for(i = 0; i < M; i++) {
                gsl_vector_set(diag_covariance, i, sqrt((*gauss->covariance)(i, i)));
            }
            gsl_vector* pca_v = gsl_vector_alloc(2);
            gsl_blas_dgemv(CblasTrans, 1.0, pca_m, diag_covariance, 0.0, pca_v);
            gsl_vector_free(diag_covariance);
            Vector* v = new Vector(3, 0);
            (*v)(0, gsl_vector_get(pca_v, 0));
            (*v)(1, gsl_vector_get(pca_v, 1));
            gsl_vector_free(pca_v);
            last_gauss_var_pca.append(v);
        }
        gsl_matrix_free(pca_m);
    } else {
        //ak niesu ziadne data
        N = selected_gaussians.size();
        if(N > 0) {
            gsl_matrix* m = gsl_matrix_alloc(M, N);
            for(git = selected_gaussians.begin(); git != selected_gaussians.end(); git++) {
                Vector* v = (*git)->mean;
                gsl_vector* vv = v->get_vector();
                gsl_matrix_set_col(m, i++, vv);
            }
            gsl_matrix* pca_m = gsl_pca(m, 2);

            for(it = last_gauss_pos_pca.begin(); it < last_gauss_pos_pca.end(); it++) {
                delete *it;
                *it = NULL;
            }
            last_gauss_pos_pca.resize(0);

            for(it = last_gauss_var_pca.begin(); it < last_gauss_var_pca.end(); it++) {
                delete *it;
                *it = NULL;
            }
            last_gauss_var_pca.resize(0);

            //vynasobi novu bazu s maticou dat
            gsl_matrix* data_pca = gsl_matrix_alloc(2, N);
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, pca_m, m, 0.0, data_pca);
            gsl_matrix_free(m);

            //vypocita vysku, sirku a priemer
            gsl_vector_view rx = gsl_matrix_row(data_pca, 0);
            pca_width = gsl_vector_max(&rx.vector) - gsl_vector_min(&rx.vector);
            gsl_vector_view ry = gsl_matrix_row(data_pca, 1);
            pca_height = gsl_vector_max(&ry.vector) - gsl_vector_min(&ry.vector);
            gsl_vector* pca_mean = gsl_vector_alloc(2);
            gsl_vector_set(pca_mean, 0, (gsl_vector_max(&rx.vector) + gsl_vector_min(&rx.vector)) / 2.0);
            gsl_vector_set(pca_mean, 1, (gsl_vector_max(&ry.vector) + gsl_vector_min(&ry.vector)) / 2.0);
            for(i = 0; i < N; i++) {
                gsl_vector_view data_col = gsl_matrix_column(data_pca, i);
                gsl_vector_sub(&data_col.vector, pca_mean);
            }
            gsl_vector_free(pca_mean);
            for(i = 0; i < N; i++) {
                Vector* v = new Vector(3, 0);
                (*v)(0, gsl_matrix_get(data_pca, 0, i));
                (*v)(1, gsl_matrix_get(data_pca, 1, i));
                last_gauss_pos_pca.append(v);
            }
            gsl_matrix_free(data_pca);

            //prepocita variancie
            for(git = selected_gaussians.begin(); git != selected_gaussians.end(); git++) {
                gsl_vector* diag_covariance = gsl_vector_alloc(M);
                Gaussian* gauss = *git;
                for(i = 0; i < M; i++) {
                    gsl_vector_set(diag_covariance, i, sqrt((*gauss->covariance)(i, i)));
                }
                gsl_vector* pca_v = gsl_vector_alloc(2);
                gsl_blas_dgemv(CblasTrans, 1.0, pca_m, diag_covariance, 0.0, pca_v);
                gsl_vector_free(diag_covariance);
                Vector* v = new Vector(3, 0);
                (*v)(0, gsl_vector_get(pca_v, 0));
                (*v)(1, gsl_vector_get(pca_v, 1));
                gsl_vector_free(pca_v);
                last_gauss_var_pca.append(v);
            }
            gsl_matrix_free(pca_m);
        } else {
            for(it = last_gauss_pos_pca.begin(); it < last_gauss_pos_pca.end(); it++) {
                delete *it;
                *it = NULL;
            }
            last_gauss_pos_pca.resize(0);

            for(it = last_gauss_var_pca.begin(); it < last_gauss_var_pca.end(); it++) {
                delete *it;
                *it = NULL;
            }
            last_gauss_var_pca.resize(0);
        }
    }
    set_wh(screen_width, screen_height);
};

void StreamArea::calc_data_gauss()
{
    double prob;
    List<Vector*>::iterator it;
    set<Gaussian*>::iterator git;
    uint i = 0;
    List<Model*>::iterator mit;
    List<State*>::iterator sit;
    List<Gaussian*>::iterator gitt;
    int index = modelset->stream_areas.index(this);
    for(mit = modelset->models.begin(); mit < modelset->models.end(); mit++) {
        for(sit = (*mit)->states.begin(); sit < (*mit)->states.end(); sit++) {
            for(gitt = (*sit)->streams[index]->gaussians.begin(); gitt < (*sit)->streams[index]->gaussians.end(); gitt++) {
                (*gitt)->my_data.resize(0);
            }
        }
    }

    for(it = data.begin(); it < data.end(); it++, i++) {
        Gaussian* gmax = NULL;
        double probmax = -DBL_MAX;
        for(git = selected_gaussians.begin(); git != selected_gaussians.end(); git++) {
            prob = (*git)->probability(*it);
            if(probmax < prob) {
                probmax = prob;
                gmax = *git;
            }
        }
        if(gmax != NULL) {
            gmax->my_data.append(i);
        }
    }
};

Vector* StreamArea::get_data(uint index)
{
    return data[index];
};

/*----------------StreamArea----------------*/

/*------------------FileData----------------*/

FileData::FileData()
{
    selected = true;
    maxprob = -DBL_MAX;
    model = NULL;
    word = "";
};

FileData::FileData(string w, List<List<Vector*> > list)
{
    selected = true;
    maxprob = -DBL_MAX;
    model = NULL;
    word = w;
    data = list;
};

List<Vector*> FileData::operator[](uint index)
{
    List<Vector*> prob_list;
    prob_list.resize(0);
    for(uint k = 0; k < data.size(); k++) {
        prob_list.append(data[k][index]);
    }
    return prob_list;
};

/*------------------FileData----------------*/

/*------------------ModelSet----------------*/

ModelSet::ModelSet(): HMMLab_Object("modelset", MODELSET), streams_size(1)
{
    stream_areas.append(new StreamArea(this, true));
};

ModelSet::ModelSet(string name): HMMLab_Object(name, MODELSET), streams_size(1)
{
    stream_areas.append(new StreamArea(this, true));
};

ModelSet::ModelSet(string filename, const char* format): HMMLab_Object("modelset", MODELSET), streams_size(1)
{
    ifstream in_stream;
    in_stream.open(filename.c_str(), fstream::in);
    load(in_stream, format);
    in_stream.close();
};

ModelSet::~ModelSet()
{
    cout << "DELETING MODELSET\n";
    uint size = models.size();
    for(uint i = 0; i < size; i++) {
        delete models[i];
    }
    size = stream_areas.size();
    for(uint i = 0; i < size; i++) {
        delete stream_areas[i];
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
                    for(uint i = 0; i < streams_size; i++) {
                        stream_areas.append(new StreamArea(this, true));
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
        uint size = streams_distribution.size();
        for(uint i = 0; i < size; i++) {
            out_stream << streams_distribution[i];
            if(i < streams_distribution.size() - 1) {
                out_stream << ' ';
            }
        }
        out_stream << endl << "<VECSIZE> " << sum(streams_distribution);
        size = vecsize_tags.size();
        for(uint i = 0; i < size; i++) {
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
        size = keys.size();
        for(uint i = 0; i < size; i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == SVECTOR && obj->ref_num > 1) {
                out_stream << "~u \"" << obj->name << '"' << endl;
                static_cast<SVector*>(obj)->save(out_stream, HTK_FORMAT);
            }
        }
        for(uint i = 0; i < size; i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == SMATRIX && obj->ref_num > 1) {
                out_stream << "~v \"" << obj->name << '"' << endl;
                static_cast<SMatrix*>(obj)->save(out_stream, HTK_FORMAT);
            }
        }
        for(uint i = 0; i < size; i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == TRANSMATRIX && obj->ref_num > 1) {
                out_stream << "~t \"" << obj->name << '"' << endl;
                uint size2 = 0;
                for(uint j = 0; j < size; j++) {
                    HMMLab_Object* om = objects_dict[keys[j]];
                    if(om->type == MODEL) {
                        Model* m = static_cast<Model*>(om);
                        if(m->trans_mat == obj) {
                            size2 = m->states.size();
                            break;
                        }
                    }
                }
                assert(size2 != 0);
                static_cast<TransMatrix*>(obj)->save(out_stream, HTK_FORMAT, size2);
            }
        }
        for(uint i = 0; i < size; i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == GAUSSIAN && obj->ref_num > 1) {
                out_stream << "~m \"" << obj->name << '"' << endl;
                static_cast<Gaussian*>(obj)->save(out_stream, HTK_FORMAT);
            }
        }
        for(uint i = 0; i < size; i++) {
            HMMLab_Object* obj = objects_dict[keys[i]];
            if(obj->type == STATE && obj->ref_num > 1) {
		    cout << "savujem\n";
                out_stream << "~s \"" << obj->name << '"' << endl;
                static_cast<State*>(obj)->save(out_stream, HTK_FORMAT);
            }
        }

        //uklada modely
        size = models.size();
        for(uint i = 0; i < size; i++) {
            Model* model = models[i];
            out_stream << "~h \"" << model->name << '"' << endl << "<BEGINHMM>" << endl;
            model->save(out_stream, HTK_FORMAT);
            out_stream << "<ENDHMM>" << endl;
        }
    } else if(!strcmp(format, XML_FORMAT)) {
    };
};

void ModelSet::create_cfg(string filename)
{
    uint d = 1, e = 1;
    ofstream cfgfile(("/dev/shm/" + filename).c_str());
    cfgfile << "SOURCEKIND = WAVEFORM" << endl;
    cfgfile << "SOURCEFORMAT = WAV" << endl;
    cfgfile << "TARGETKIND = ";
    if(vecsize_tags.index(mfcc_e_d_a) != -1) {
        cfgfile << "MFCC_E_D_A";
        d = 3;
    } else if(vecsize_tags.index(mfcc_e_d) != -1) {
        cfgfile << "MFCC_E_D";
        d = 2;
    } else if(vecsize_tags.index(mfcc_e_d) != -1) {
        cfgfile << "MFCC_E";
    } else {
        cfgfile << "MFCC";
        e = 0;
    }
    cfgfile << endl << "TARGETRATE = 100000" << endl;
    cfgfile << "TARGETFORMAT = HTK" << endl;
    cfgfile << "WINDOWSIZE = 250000" << endl;
    cfgfile << "ZMEANSOURCE = T" << endl;
    cfgfile << "USEHAMMING = T" << endl;
    cfgfile << "PREEMCOEF = 0.97" << endl;
    cfgfile << "ADDDITHER = 1.1" << endl;
    cfgfile << "NUMCHANS = 40" << endl;
    cfgfile << "CEPLIFTER = 22" << endl;
    cfgfile << "ENORMALISE = F" << endl;
    cfgfile << "NUMCEPS = " << (dimension / d) - e << endl;
    cfgfile.close();
};

void ModelSet::load_data(uint length, string* filenames)
{
    create_cfg("hmmlab.cfg");
    for(uint i = 0; i < length; i++) {
        stringstream cmd(stringstream::in | stringstream::out);
        cmd << "HList -C /dev/shm/hmmlab.cfg -r -n " << streams_size << " " << filenames[i];
        stringstream str_data(stringstream::in | stringstream::out);
        str_data << execute(cmd.str());
        List<Vector*> list_data;
        double num;
        while(str_data.good()) {
            Vector* vec = new Vector(dimension);
            uint i;
            for(i = 0; i < dimension; i++) {
                str_data >> scientific >> num;
                (*vec)(i, num);
            }
            list_data.append(vec);
        }
        add_data(filenames[i], list_data);
    }
    unlink("/dev/shm/hmmlab.cfg");
    reset_pos_gauss();
};

void ModelSet::add_data(string filename, List<Vector*> d)
{
    uint start = 0, size;
    for(uint i = 0; i < streams_size; i++) {
        List<Vector*> list;
        size = d.size();
        for(uint k = 0; k < size; k++) {
            Vector* vec = new Vector(streams_distribution[i], 0);
            for(uint l = 0; l < streams_distribution[i]; l++) {
                (*vec)(l, (*d[k])[start + l]);
            }
            list.append(vec);
        }
        start += streams_distribution[i];

        stream_areas[i]->add_data(list);
        if(files_data.find(filename) == files_data.end()) {
            List<List<Vector*> > fdata;
            for(uint k = 0; k < streams_size; k++) {
                List<Vector*> vlist;
                fdata.append(vlist);
            }
            ifstream in_stream;
            string str = filename.substr(0, filename.length() - 3) + "lab";
            string word = "";
            if(access(str.c_str(), F_OK) != -1) {
                in_stream.open(str.c_str(), fstream::in);
                in_stream >> word;
                while(!isallalpha(word)) {
                    in_stream >> word;
                }
            } else {
                uint found = filename.find_last_of('/');
                string tmp = filename.substr(found + 1, filename.length() - found);
                for(uint j = 0; j < tmp.length(); j++) {
                    if(isalpha(tmp[j])) {
                        word += tmp[j];
                    } else {
                        break;
                    }

                }
            }
            files_data[filename] = new FileData(word, fdata);
        }
        files_data[filename]->data[i] = list;
    }
};

void ModelSet::reset_pos_gauss()
{
    List<StreamArea*>::iterator it;

    for(it = stream_areas.begin(); it < stream_areas.end(); it++) {
        (*it)->reset_pos_gauss();
    }
    uint i = 0, size;
    set<Model*>::iterator mit;
    size = drawarea_models.size();
    Dict<string, FileData*>::iterator fdit;
    for(fdit = files_data.begin(); fdit != files_data.end(); fdit++) {
        fdit->second->model = NULL;
        fdit->second->maxprob = -DBL_MAX;
    }
    Model** array = new Model*[size];
    for(mit = drawarea_models.begin(); mit != drawarea_models.end(); mit++, i++) {
        array[i] = *mit;
    }
    run_threads<Model*>(size, array, run_model_viterbi);
    delete[] array;
};

bool ModelSet::is_selected(Model* m, int index)
{
    State* s = m->states[index];
    bool sel = true;
    List<Stream*>::iterator it;
    List<Gaussian*>::iterator git;
    List<StreamArea*>::iterator sit;
    for(it = s->streams.begin(); it < s->streams.end(); it++) {
        for(git = (*it)->gaussians.begin(); git < (*it)->gaussians.end(); git++) {
            bool in = false;
            for(sit = stream_areas.begin(); sit < stream_areas.end(); sit++) {
                in |= (*sit)->selected_gaussians.count(*git) != 0;
            }
            sel &= in;
        }
    }
    return sel;
};

bool ModelSet::is_selected(Gaussian* g)
{
    List<StreamArea*>::iterator sit;

    for(sit = stream_areas.begin(); sit < stream_areas.end(); sit++) {
        if((*sit)->selected_gaussians.count(g) != 0) {
            return true;
        }
    }
    return false;
};

uint ModelSet::selected_gaussians_count()
{
    uint count = 0;
    List<StreamArea*>::iterator sit;
    for(sit = stream_areas.begin(); sit < stream_areas.end(); sit++) {
        count += (*sit)->selected_gaussians.size();
    }
    return count;
};

uint ModelSet::loaded_data_count()
{
    uint count = 0;
    List<StreamArea*>::iterator sit;
    for(sit = stream_areas.begin(); sit < stream_areas.end(); sit++) {
        count += (*sit)->data.size();
    }
    return count;
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

void ModelSet::gnuplot_2D(uint stream_index, uint dim)
{
    assert(stream_index >= 0 && stream_index < streams_size);
    uint j = 0;
    Gaussian* gauss;
    set<Gaussian*>::iterator it;
    gnuplot_ctrl* h = gnuplot_init();
    char buffer[256];
    gnuplot_cmd(h, "set style fill transparent solid 0.50 noborder");
    sprintf(buffer, "set title \"Gaussian mixture dimension %d stream %d\"", dim, stream_index);
    gnuplot_cmd(h, buffer);
    gnuplot_cmd(h, "set key inside left top vertical Left reverse enhanced autotitles nobox");
    gnuplot_cmd(h, "set key noinvert samplen 1 spacing 1 width 0 height 0 ");
    gnuplot_cmd(h, "set style function filledcurves y1=0");
    gnuplot_cmd(h, "set tmargin 2");
    gnuplot_cmd(h, "Gauss(x, mu, sigma) =  1./(sqrt(2*pi*sigma)) * exp( -(x-mu)**2 / (2*sigma) )");
    gnuplot_cmd(h, buffer);
    gnuplot_cmd(h, "unset key");
    double max = -DBL_MAX;
    double min = DBL_MAX;
    j = 0;
    stringstream cmd(stringstream::in | stringstream::out);
    stringstream cmd_sum(stringstream::in | stringstream::out);
    set<Gaussian*>* sel_gauss = &stream_areas[stream_index]->selected_gaussians;
    for(it = sel_gauss->begin(); it != sel_gauss->end(); it++) {
        gauss = *it;
        double mu = (*(gauss->mean))[dim];
        double sigma = (*(gauss->covariance))(dim, dim);
        double sigma_sqrt = sqrt(sigma);
        if(mu + sigma_sqrt * 2  > max) {
            max = mu + sigma_sqrt * 2;
        }
        if(mu - sigma_sqrt * 2 < min) {
            min = mu - sigma_sqrt * 2;
        }
        cmd << "Gauss(x," << scientific << mu << ',' << scientific << sigma << "), ";
        cmd_sum << "Gauss(x," << scientific << mu << ',' << scientific << sigma << ")";
        if(++j < sel_gauss->size()) {
            cmd_sum << " + ";
        }
    }
    cmd << cmd_sum.str();

    if(stream_areas[stream_index]->data.size() > 0) {
        stream_areas[stream_index]->save_data_pos_2D(dim, "data");
        cmd << ", \"/dev/shm/data\"";
    }

    stringstream ccmd(stringstream::in | stringstream::out);
    ccmd << "plot [" << scientific << min << ':' << scientific << max << "] " << cmd.str();
    string str = ccmd.str();
    char* writable = new char[str.size() + 1];
    copy(str.begin(), str.end(), writable);
    writable[str.size()] = '\0';
    gnuplot_cmd(h, writable);
    delete[] writable;
    //gnuplot_close(h);
};

string ModelSet::get_unique_name(string prefix)
{
    char* buffer = new char[prefix.size() + 100];
    bool got_name = false;
    uint i;
    for(i = 0; !got_name; i++) {
        sprintf(buffer, "%s%d", prefix.c_str(), i);
        got_name = objects_dict.count(buffer) <= 0;
    }
    string ret = buffer;
    delete buffer;
    return ret;
};

List<Model*> ModelSet::get_models_with_gaussian(Gaussian* g)
{
    List<Model*> ret;
    List<Model*>::iterator mit;
    List<State*>::iterator sit;
    List<Stream*>::iterator strit;
    for(mit = models.begin(); mit < models.end(); mit++) {
        bool has_gauss = false;
        for(sit = (*mit)->states.begin(); sit < (*mit)->states.end(); sit++) {
            for(strit = (*sit)->streams.begin(); strit < (*sit)->streams.end(); strit++) {
                has_gauss = (*strit)->gaussians.index(g) != -1;
                if(has_gauss) {
                    break;
                }
            }
            if(has_gauss) {
                break;
            }
        }
        if(has_gauss) {
            ret.append(*mit);
        }
    }
    return ret;
};

bool ModelSet::gauss_cluster(List<Gaussian*> gaussians, List<Vector*> data)
{
    double Psum, new_log_likelihood = 0.0, log_likelihood = -DBL_MAX;
    uint i, j;
    uint dimension = streams_distribution[gaussians[0]->index_distribution];
    uint k = gaussians.size();
    uint n = data.size();
    gsl_vector_view col;
    gsl_vector* d = gsl_vector_alloc(dimension);
    gsl_vector* tmp;
    gsl_matrix* mat, *tmpmat;

    gsl_vector* pi = gsl_vector_alloc(k);
    gsl_vector_set_all(pi, 1.0 / k);
    gsl_matrix* prob = gsl_matrix_alloc(n, k);
    gsl_vector* f = gsl_vector_alloc(n);
    gsl_matrix* P = gsl_matrix_alloc(n, k);

    while(log_likelihood < new_log_likelihood) {
        log_likelihood = new_log_likelihood;
        //pravdepodobnosti, kazdy s kazdym
        for(i = 0; i < n; i++) {
            for(j = 0; j < k; j++) {
                gsl_matrix_set(prob, i, j, hmmlab_exp(gaussians[j]->probability(data[i])));
            }
        }

        //vypocita f
        gsl_blas_dgemv(CblasNoTrans, 1.0, prob, pi, 0.0, f);

        //vypocita P
        gsl_matrix_memcpy(P, prob);
        for(j = 0; j < n; j++) {
            col = gsl_matrix_row(P, j);
            gsl_vector_mul(&col.vector, pi);
        }
        for(i = 0; i < k; i++) {
            col = gsl_matrix_column(P, i);
            gsl_vector_div(&col.vector, f);
        }

        //vypocita pi
        for(i = 0; i < k; i++) {
            col = gsl_matrix_column(P, i);
            gsl_vector_set(pi, i, gsl_stats_mean(col.vector.data, col.vector.stride, col.vector.size));
        }

        //prepocita stredy gaussianov
        for(i = 0; i < k; i++) {
            Psum = 0.0;
            tmp = gaussians[i]->mean->get_vector();
            gsl_vector_set_all(tmp, 0.0);
            for(j = 0; j < n; j++) {
                Psum += gsl_matrix_get(P, j, i);
                gsl_vector_memcpy(d, data[j]->get_vector());
                gsl_vector_scale(d, gsl_matrix_get(P, j, i));
                gsl_vector_add(tmp, d);
            }
            gsl_vector_scale(tmp, 1.0 / Psum);
        }

        //prepocita kovariancne matice
        for(i = 0; i < k; i++) {
            Psum = 0.0;
            tmp = gaussians[i]->mean->get_vector();
            mat = gaussians[i]->covariance->get_matrix();
            gsl_matrix_set_all(mat, 0);
            for(j = 0; j < n; j++) {
                Psum += gsl_matrix_get(P, j, i);
                gsl_vector_memcpy(d, data[j]->get_vector());
                gsl_vector_sub(d, tmp);
                gsl_blas_dger(gsl_matrix_get(P, j, i), d, d, mat);
            }
            gsl_matrix_scale(mat, 1.0 / Psum);
        }

        //diagonalizuje
        for(i = 0; i < k; i++) {
            mat = gaussians[i]->covariance->get_matrix();
            gsl_vector_view diag = gsl_matrix_diagonal(mat);
            gsl_vector_memcpy(d, &diag.vector);
            for(j = 0; j < dimension; j++) {
                if(gsl_vector_get(d, i) < COVMIN) {
                    gsl_vector_set(d, i, COVMIN);
                }
            }
            gsl_matrix_set_all(mat, 0.0);
            gsl_vector_memcpy(&diag.vector, d);
            tmpmat = gaussians[i]->inv_covariance->get_matrix();
            for(j = 0; j < dimension; j++) {
                gsl_matrix_set(tmpmat, j, j, 1.0 / gsl_matrix_get(mat, j , j));
            }
            gaussians[i]->calc_gconst();
        }
        new_log_likelihood = 0;
        for(i = 0; i < n; i++) {
            new_log_likelihood += hmmlab_log(gsl_vector_get(f, i));
        }
    }

    gsl_vector_free(d);
    gsl_vector_free(pi);
    gsl_vector_free(f);
    gsl_matrix_free(prob);
    gsl_matrix_free(P);
    List<StreamArea*>::iterator it;
    for(it = stream_areas.begin(); it < stream_areas.end(); it++) {
        (*it)->calc_data_gauss();
    }
    reset_pos_gauss();

    return true;
};

bool ModelSet::gauss_push(bool out, Gaussian* g1, Gaussian* g2)
{
    Vector dist1 = (*g1->mean - *g2->mean) * (out ? GAUSS_PUSH : - GAUSS_PUSH);
    Vector dist2 = (*g2->mean - *g1->mean) * (out ? GAUSS_PUSH : - GAUSS_PUSH);
    SVector* m1 = new SVector(g1->mean->name, this, g1->mean->size(), 0.0);
    SVector* m2 = new SVector(g2->mean->name, this, g2->mean->size(), 0.0);
    *m1 = *g1->mean;
    *m2 = *g2->mean;
    gsl_vector_add(m1->get_vector(), dist1.get_vector());
    gsl_vector_add(m2->get_vector(), dist2.get_vector());
    g1->mean->dec_ref_num();
    g2->mean->dec_ref_num();
    g1->mean = m1;
    g2->mean = m2;
    List<StreamArea*>::iterator it;
    for(it = stream_areas.begin(); it < stream_areas.end(); it++) {
        (*it)->calc_data_gauss();
    }
    reset_pos_gauss();
    return true;
};

void ModelSet::drawarea_models_append(Model* m)
{
    drawarea_models.insert(m);
    uint i = 0;
    set<Model*>::iterator mit;
    Dict<string, FileData*>::iterator fdit;
    for(fdit = files_data.begin(); fdit != files_data.end(); fdit++) {
        fdit->second->model = NULL;
        fdit->second->maxprob = -DBL_MAX;
    }
    Model** array = new Model*[drawarea_models.size()];
    for(mit = drawarea_models.begin(); mit != drawarea_models.end(); mit++, i++) {
        array[i] = *mit;
    }
    run_threads<Model*>(drawarea_models.size(), array, run_model_viterbi);
    delete[] array;
};

void ModelSet::select_data(string filename)
{
    uint i, j, size, size2;
    Vector* vec;
    files_data[filename]->selected = true;
    List<List<Vector*> > vectors = files_data[filename]->data;
    size = vectors.size();
    for(i = 0; i < size; i++) {
        size2 = vectors[i].size();
        for(j = 0; j < size2; j++) {
            vec = vectors[i][j];
            if(stream_areas[i]->data.index(vec) == -1) {
                stream_areas[i]->data.append(vec);
            }
        }
        stream_areas[i]->calc_edge_len();
        stream_areas[i]->calc_data_gauss();
    }
    reset_pos_gauss();
};

void ModelSet::unselect_data(string filename)
{
    uint i, j, size, size2;
    files_data[filename]->selected = false;
    List<List<Vector*> > vectors = files_data[filename]->data;
    size = vectors.size();
    for(i = 0; i < size; i++) {
        size2 = vectors[i].size();
        for(j = 0; j < size2; j++) {
            stream_areas[i]->data.remove_value(vectors[i][j]);
        }
        stream_areas[i]->calc_edge_len();
        stream_areas[i]->calc_data_gauss();
    }
    reset_pos_gauss();
};

void ModelSet::train_model(Model* m)
{
    List<FileData*> data;
    Dict<string, FileData*>::iterator it;
    for(it = files_data.begin(); it != files_data.end(); it++) {
        if(it->second->selected) {
            data.append(it->second);
        }
    }
    if(!data.empty()) {
        m->train(data);
    }
};
/*------------------ModelSet----------------*/

#endif
