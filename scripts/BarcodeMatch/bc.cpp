#include "bc.h"
#include "edit.h"
#include "normal.h"
using namespace std;


vector<int> kmer_include(string seq, set<string> dic) {
    vector<int> sc;
    int id = 0;
    for (set<string>::iterator it = dic.begin(); it != dic.end(); it++) {
        if (seq.find(*it) != seq.npos) {
            sc.push_back(id);
        }
        id++;
    }
    return(sc);
}

vector<vector<int> > barcodes_cos_vector(vector<string> barcodes, set<string> dic) {
    vector<vector<int> > bar_kmer_vec;
    int n = barcodes.size();

    for (int i = 0; i < n; i++) {
        int id = 0;
        vector<int> temp_vec = kmer_include(barcodes[i],dic);
        bar_kmer_vec.push_back(temp_vec);
    }
    return(bar_kmer_vec);
}

double cos_sim(vector<int> a, vector<int> b) {
    if (a.size() == 0 || b.size() == 0) {
        return(0);
    }
    vector<int> inter;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(inter));

    return(inter.size() / (sqrt(a.size()) * sqrt(b.size())));
}

bool cos_sim_comp(pair <int, double> cos1, pair <int, double> cos2) {
    return(cos1.second > cos2.second);
}

vector<string> barcode_cand_cos(string seq, vector<string> barcodes, set<string> dic, vector<vector<int> > index,
    int top, double thresh) {

    int bar_size = barcodes.size();

    vector<pair<int, double> > bar_cand;

    vector<int> sc = kmer_include(seq, dic);
    for (int j = 0; j < bar_size; j++) {
        double cos = cos_sim(sc,index[j]);
        if (cos >= thresh) {
            bar_cand.push_back(pair<int, double>(j, cos));
        }
    }
    sort(bar_cand.begin(), bar_cand.end(), cos_sim_comp);
    int count = bar_cand.size();
    int limit = min(top, count);

    vector<string> out;

    for (int i = 0; i < limit; i++) {
        out.push_back(barcodes[bar_cand[i].first]);
    }
    return(out);
}

vector<int> pos_filter(vector<double> start, vector<int> edit) {
    int num = start.size();

    vector<double> correct;

    int edit_count[16] = {}, offset_count[30] = {};
    int max_edit = 0, max_offset = 0;

    double edit_sum = 0;

    for (int i = 0; i < num; i++) {
        if (edit[i] == 0) {
            correct.push_back(start[i]);
        }
        if (edit[i] > max_edit) {
            max_edit = edit[i];
        }
        edit_count[edit[i]]++;
        edit_sum += edit[i];
    }

    double edit_mean = edit_sum / num;

    double mu = mean(correct);

    for (int i = 0; i < num; i++) {
        start[i] = round(abs(start[i] - mu));
        if (start[i] > max_offset) {
            max_offset = start[i];
        }
        offset_count[int(start[i])]++;
    }

    double edit_ratio[16] = {}, offset_ratio[30] = {};
    for (int i = max_edit; i >= 0; i--) {
        double count = 0;
        for (int j = i; j <= max_edit; j++) {
            count += edit_count[j];
        }
        edit_ratio[i] = count / double(num);
    }
    for (int i = max_offset; i >= 0; i--) {
        double count = 0;
        for (int j = i; j <= max_offset; j++) {
            count += offset_count[j];
        }
        offset_ratio[i] = count / double(num);
    }

    vector<pair<int, int> > preserve;
    for (int i = 0; i <= max_edit; i++) {
        for (int j = 0; j <= max_offset; j++) {
            if (pow(edit_ratio[i], edit_mean) * pow(offset_ratio[j], 0.5) >= 0.05 || i == 0) {
                preserve.push_back(pair<int, int>(i, j));
            }
        }
    }

    vector<int> preserve_id;
    for (int i = 0; i < preserve.size(); i++) {
        for (int j = 0; j < num; j++) {
            if (edit[j] == preserve[i].first && start[j] == preserve[i].second) {
                preserve_id.push_back(j);
            }
        }
    }
    sort(preserve_id.begin(), preserve_id.end());
    return(preserve_id);
}

void barcodeMatch(vector<string> seq, vector<string> barcodes, vector<string> readnames,
    double mu, double sigma, double sigma_start, int k, int batch,
    int top, double cos_thresh, double alpha, int edit_thresh,
    string outname) {

    printf("start to build index\n");
    set<string> dic = kmer(barcodes, k);
    vector<vector<int> > index = barcodes_cos_vector(barcodes, dic);
    printf("index finished!\n");

    int seq_size = seq.size();
    int bar_size = barcodes.size();

    int bar_len = barcodes[1].length();

    int times = seq_size / batch;
    int seq_s, seq_e;

    ofstream OutFile(outname);
    vector<int> result_id;
    vector<string> result_bar;
    vector<double> result_s;
    vector<int> result_edit;

    for (int i = 1; i <= times + 1; i++) {
        seq_s = batch * (i - 1);
        seq_e = batch * i;

        if (i == times + 1) {
            if (seq_size % batch) {
                seq_e = seq_size;
            }
            else {
                break;
            }
        }
        double interval_s = qnorm(alpha / 2, mu, sigma) < 0 ? 0 : qnorm(alpha / 2, mu, sigma);
        double interval_e = qnorm(1 - alpha / 2, mu, sigma) + bar_len;

        for (int j = seq_s; j < seq_e; j++) {
            if (seq[j].length() < interval_s+bar_len-1) {
                continue;
            }
            string sub_seq = seq[j].substr(interval_s, interval_e - interval_s + 1);

            if (sub_seq.length() < bar_len) {
                continue;
            }
            else {
                vector<string> retain_barcodes = barcode_cand_cos(sub_seq, barcodes, dic, index, top, cos_thresh);
                int retain_barcodes_size = retain_barcodes.size();

                if (retain_barcodes_size == 0) {
                    continue;
                }
                int now = bar_len;
                pair<int, int> best(now, -100);
                int bar_id = 0;
                bool flag = 0;
                for (int p = 0; p < retain_barcodes_size; p++) {
                    pair<int, int> temp = minEditDist(sub_seq, retain_barcodes[p]);
                    if (temp.first < now) {
                        flag = 1;
                        bar_id = p;
                        best = temp;
                        now = temp.first;
                    }
                    else if (temp.first == now) {
                        if (abs(temp.second - mu) < abs(best.second - mu)) {
                            best = temp;
                            bar_id = p;
                        }
                    }
                }

                if (!flag) {
                    continue;
                }
                else if (best.first <= edit_thresh) {
                    result_id.push_back(j);
                    result_bar.push_back(retain_barcodes[bar_id]);
                    result_s.push_back(interval_s + best.second);
                    result_edit.push_back(best.first);
                }
            }

        }
        mu = mean(result_s);
        sigma = update_sigma(result_s, sigma_start);
    }

    
    vector<int> preserve_id = pos_filter(result_s, result_edit);
    int preserve_num = preserve_id.size();
    cout << "There are " << preserve_num << " sequences identified with a barcode" << endl;

    for (int i = 0; i < preserve_num; i++) {
        OutFile << result_id[preserve_id[i]] << " " << readnames[result_id[preserve_id[i]]] << " "
            << result_bar[preserve_id[i]] << " " << result_s[preserve_id[i]] << " " << result_edit[preserve_id[i]] << endl;
    }
    
    OutFile.close();
}