#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
using namespace std;
// author: Jing
// time: April 20 2019
// run: clang++ -std=c++14 read_matrix.cpp -o rm && ./rm
// 1-s2.0-S1074761318301213-mmc2.xlsx TCGA study
//mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriverLoose_aggregated_matrix.tsv
// put files in ../DATA/

inline bool exists_test0 (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}
bool isFloat( string myString ) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}
// read the mc3 file and find index of one patient
vector<int> removeDuplicate(){
   string filename("../mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriverLoose_aggregated_matrix.tsv");
  
  ifstream datafile;
  datafile.open(filename);

  ofstream myfile("matrix.txt");
  
  string fn, tmp;
  getline(datafile, fn);
  datafile.close();
  
  stringstream ss(fn);
  unordered_map<string, vector<int>> dict;
  int i = 0;
  while (ss){
     ss >> tmp;
     dict[tmp.substr(0,12)].emplace_back(i);
     ++i;
  }
  ofstream wfile("pid_index.txt"); // remove Hugo_Symbol
  vector<int> id;
  for( auto iter = dict.begin(); iter != dict.end(); ++iter ){
     cout << iter->first << '\t' << iter->second[0] << endl;
     wfile << iter->first << '\t' << iter->second[0] << endl;
     id.emplace_back(iter->second[0]);
  }
  wfile.close();
  return id;
}
// pid_index provides the column of mc3 that we need
//void readFile(){
//  string tmp, tmp_ss;
//  int i (0);
//  ifstream pid_file;
//  pid_file.open( "pid_index.txt" );
//  unordered_set<string> patient;
//  set<int> index;
//  while( getline(pid_file, tmp) ){
//     stringstream ss( tmp );
//     i = 0;
//     while( ss >> tmp_ss ){
//       if(i == 0) {
//            patient.emplace( tmp_ss );
//       }else{
//            index.emplace( stoi(tmp_ss) );
//            if( stoi(tmp_ss) < 1)
//               cout << tmp << " error" << endl;
//       }
//       ++i;
//     }
//     ss.str( "" );
//  }
//  pid_file.close();
//  cout << "patient : " << patient.size() << '\t' << "index: " << index.size() << endl;
//  cout << *index.begin() << '\t' << *index.rbegin() << endl;
//  ifstream datafile;
//  string filename("../DATA/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriverLoose_aggregated_matrix.tsv");
//  datafile.open(filename); 
// // ofstream myfile( "matrix_gene_pid.txt" ); //"matrix_gene_pid.txt"
//  vector<vector<string>> matrix_pid_fea (index.size(), vector<string>( 300, " "));
//   cout << matrix_pid_fea.size() << '\t' << matrix_pid_fea[0].size() << endl;
//  int fea_ind (0), samp_ind(0), new_row(0); // k is column of mc3, j is feature
//  cout << "hello" << endl;
//  while( getline(datafile, tmp)){
//    stringstream ss(tmp);
//    samp_ind = 0;
//    new_row = 0;
//    //cout << tmp << endl;
//    while( ss ){ // i is the index of column, j is row feature index    
//        ss >> tmp_ss;
//        string trans_string;
//        if ( index.count(samp_ind) > 0 ){
//    //     cout << tmp_ss << '\t';
//           if (tmp_ss.size() > 5 && tmp_ss.substr(0,4) == "TCGA"){
////             myfile << tmp_ss.substr(0,12);
//               trans_string = tmp_ss.substr(0, 12);
//           }else if (tmp_ss.size() >= 2 && tmp_ss.substr(0,2) == "NA"){
////             myfile << 1;
//               trans_string = "1";
//           }else if (tmp_ss.size() >= 2 && tmp_ss.substr(0,2) == "p."){
////             myfile << 0;
//              // cout << tmp_ss << '\t';
//               trans_string = "0";
//           }else{
////             myfile << tmp_ss;
//               trans_string = tmp_ss;
//           }
////         myfile << '\t';
//         //  cout << samp_ind << '\t' << fea_ind<<'\t';
//           matrix_pid_fea[new_row++][fea_ind] = trans_string;
//       }
//      // cout << endl;
//       ++samp_ind;
//    }
//    ++fea_ind;
//    myfile << endl;
//  }
//  myfile.close();
//  
//  ofstream transfile;
//  transfile.open( "matrix_pid_gene.txt" );
//  for(int t = 0; t < index.size(); ++t){
//     for(string s : matrix_pid_fea[t]){
//         transfile << s << '\t';
//     }
//     transfile << endl;
//  }
//  transfile.close();
//}
void read_first_col(string txt_name, string fc_name){ // the first column of mc3...tsv is the genes
  string tmp, tmp_ss;
  
  ifstream data_file( txt_name );
  ofstream gen_file( fc_name );
  int i (0);
  int j (0);
  while( getline(data_file, tmp) ){
    stringstream ss(tmp); 
    i = 0;
//    while( getline(ss, tmp_ss, ',') && i == 0){
      while(ss >> tmp_ss && i++ == 0){
       gen_file << tmp_ss << '\t' << j++;
    }
    gen_file << endl;
  }
  gen_file.close();
  data_file.close();
  cout << "Read the first column of " << txt_name << " to " << fc_name << " done \n";
}
void read_first_row(string txt_name, string fr_name){ 
  // the first row of mc3...tsv is the genes
  string tmp, tmp_ss;
  
  ifstream data_file( txt_name );
  ofstream gen_file( fr_name );
  int i (0);
  int j (0);
  while( getline(data_file, tmp) && i == 0){
    stringstream ss(tmp); 
//    while( getline(ss, tmp_ss, ',') && i == 0){
    while( getline(ss, tmp_ss, '\t')){
       gen_file << tmp_ss << '\t' << j++ << '\n';
    }
    ++i;
  }
  gen_file.close();
  data_file.close();
  cout << "Read the first row of " << txt_name << " to " << fr_name << " done \n";
}
void read_col_index(string txt_name, string ft_name, int index){ 
  // the first row of mc3...tsv is the genes
  string tmp, tmp_ss;
  
  ifstream data_file( txt_name );
  if (!data_file) std::cerr << "Cound not open the" << txt_name << "file! \n";

  int i (0), j(0), k(0);
  map<string, int> dict_col;
  while( getline(data_file, tmp) ){
    stringstream ss(tmp); 
    j = 0;
    if(i>0){
        while( getline(ss, tmp_ss, '\t')){
        if(j == index){
            if(dict_col.count(tmp_ss) <=0) dict_col[tmp_ss] = k++;
        }     
        ++j;
        }
    }
    ++i;
  }
  data_file.close();
  ofstream gen_file( ft_name ); //ios::out | ios::in 
  
  for( auto iter = dict_col.begin(); iter!=dict_col.end(); ++iter){
    gen_file << iter->first << '\t' << iter->second << '\n';
  }
  gen_file.close();
  cout << "\n Summary the " << index << " column of " << txt_name << " to " << ft_name << " done \n";
}

void readTXT_symbols(string txt_name, string summ_name){
    ifstream txt_file( txt_name );
    if (!txt_file) std::cerr << "Cound not open the TXT file! \n";
    
    map<string, int> unique_symb;
    
    string tmp;
    string tmp_ss;
    while( getline(txt_file, tmp) ){
        stringstream ss(tmp);
       // while( getline(ss, tmp_ss, ',') ){
       while( ss >> tmp_ss){
            //if(isalpha(tmp_ss[0])){
                unique_symb[tmp_ss]++;
           // }
        }
    }
    txt_file.close();
    
    ofstream summ_file( summ_name );
    for(auto iter = unique_symb.begin(); iter != unique_symb.end(); ++iter){
        summ_file << iter->first << '\t' << iter->second << '\n';
    }
    summ_file.close();
    cout << "read unique symboles from " << txt_name << " to " << summ_name<< " done! \n";
}
void read_TCGA(string tcga_name, string tcga_patr, string bar_name, string ind_name, string summ_name ){
   ifstream tcga_file ( tcga_name );
   if(!tcga_file) std::cerr << "Could not open the " << tcga_name << " file! \n";
   string tmp, tmp_code;
   map<string, vector<pair<string,int>>> dict_tcga;
   int i(0);
   while( getline(tcga_file, tmp) ){
        if( i > 0 ){
        stringstream ss(tmp);
            getline(ss, tmp_code, '\t');
            dict_tcga[tmp.substr(0,12)].emplace_back( make_pair(tmp_code,i) );
        }
        ++i;
   }
   tcga_file.close();
   
   map<int, vector<string>> dict_summ;
   ofstream code_file( tcga_patr );
   ofstream bar_file( bar_name );
   ofstream ind_file( ind_name );
   for( auto iter = dict_tcga.begin(); iter != dict_tcga.end(); ++iter ){
        code_file << iter->first << '\t';
        bar_file << iter->first << '\t';
        ind_file << iter->first << '\t';
        dict_summ[ iter->second.size() ].emplace_back( iter->first );
        for( auto & code : iter->second ){
            code_file << code.first << '\t' << code.second << '\t';
            //   getline(tcgac, tmp_code, '\t');
            bar_file << code.first << '\t';
            ind_file << code.second << '\t';
        }
        code_file << endl;
        bar_file << endl;
        ind_file << endl;
   }
   code_file.close();
   bar_file.close();
   ind_file.close();
   ofstream summ_code_file( summ_name );
   for( auto code_iter = dict_summ.begin(); code_iter != dict_summ.end(); ++code_iter ){
        summ_code_file << code_iter -> first << '\t';
        for(string tmp_str : code_iter->second){
            summ_code_file << tmp_str << '\t';
        } 
        summ_code_file << endl;
        cout << code_iter->first << '\t' << code_iter->second.size() << endl;
   }
   summ_code_file.close();
   cout << "Write " << tcga_name << " to " << tcga_patr << endl;
}

void readTSV_toTXT(string tsv_name, string txt_name){
  ifstream tsv_file( tsv_name );
  if (!tsv_file) std::cerr << "Could not open the TSV file!" << std::endl;

  ofstream txt_file( txt_name );
  string tmp;
  while( getline(tsv_file, tmp) ){
     txt_file << tmp << endl;
  }
  tsv_file.close();
  txt_file.close();
  cout << "write done!" << endl;
}
void transpose_txt_file(string txt_name, string transp_name){
    ifstream txt_file(  txt_name    );
    if(!txt_file) std::cerr << "Could not open " << txt_name << " file! \n";
    
    vector<vector<string>> row_col;
    string tmp, tmp_ss;
    vector<string> vec_str;
    while(getline(txt_file, tmp)){
        stringstream tmp_ss(tmp);
        vec_str = {};
        while(tmp_ss >> tmp){
            vec_str.emplace_back(tmp);
        }
        row_col.emplace_back( vec_str );
    }
    int num_row = row_col.size();
    int num_col = row_col[0].size();
    cout << "row: " << num_row << " col: " << num_col << endl;
    ofstream transp_file(   transp_name );
    for(int i = 0; i < num_col; ++i){
        for(int j = 0; j < num_row; ++j){
            if( j != num_row-1){
                transp_file << row_col[j][i] << '\t';
            }else{
                transp_file << row_col[j][i];
            }
        }
        transp_file << '\n';
    }
    cout << "Transpose " << txt_name << " to " << transp_name << " done! \n"; 
}
void intersection_first_col (string mf_name, string mfc_name, string mt_name){ // mf_name contains the first col in mfc_name
    set<string> mf_id;
    ifstream read_file( mf_name );
    string tmp, tmp_ss;
    int i(0);
    while( getline(read_file, tmp) ){
        if( i > 0 ){
            stringstream ss(tmp);
            getline(ss, tmp_ss, '\t');
            mf_id.emplace(tmp_ss);
        }
        ++i;
    }
    cout << mf_id.size() << endl;
    read_file.close();
    ifstream mfc_file( mfc_name );
    vector<string> gene_data;
    i = 0;
    while( getline(mfc_file, tmp) ){
        if( i++ > 0 ){
            string tcga = tmp.substr(0, 12);
            if( mf_id.count(tcga) > 0 ){
                gene_data.emplace_back( tmp );     
            }         
        }else{
            gene_data.emplace_back( tmp );
        }
    }
    mfc_file.close();
    cout << gene_data.size() << endl;
    ofstream wt_file( mt_name );
    for(auto iter = gene_data.begin(); iter != gene_data.end(); ++iter){
        wt_file << *iter << '\n';
       // cout << *iter << endl;
    }
    wt_file.close();
    cout << "The intersection of " << mfc_name << " and " << mf_name << " to " << mt_name << endl;
}
void merge_two (string mf1_name, string mf2_name, string mt_name, set<int> index) {
// merge mf1_name and the columns of mf2_name  
    ifstream mf1_file( mf1_name );
    ifstream mf2_file( mf2_name );
    if(!mf1_file) std::cerr << "Could not open " << mf1_name << " file! \n";
    if(!mf2_file) std::cerr << "Could not open " << mf2_name << " file!\n";
    string tmp, tmp_ss, name;
    map<string, vector<string>> mf1_dict;
    map<int,vector<int>>summ;
    int i( 0 );
    while( getline(mf1_file, tmp) ){
        if( i > 0){
            tmp_ss = tmp.substr(0, 12);
            
            mf1_dict[tmp_ss].emplace_back(tmp);
//            cout << tmp<<endl;
        }else{
            mf1_dict["0"].emplace_back( tmp ); 
        }
        ++i;                
    }   
//    cout << mf1_dict.size() << endl;
    int j(0);
    int s(0); 
    string ss_summ;
    while( getline(mf2_file, tmp) ){
//    cout << "j: "<<j << endl;
       if ( j > 0){
        name = tmp.substr(0, 12);
//        cout << name << endl;
        if( mf1_dict.count(name) >0 ){
          stringstream ss(tmp);
          i = 0;
          while( getline(ss, tmp_ss, '\t') ){
            s = 0;
            if( index.count(i) ){ 
//            cout << "i "<<i << endl;
                for(int k = 0; k < mf1_dict[name].size(); ++k){
//                    stringstream summ(mf1_dict[name][k]);
//                    while(getline(summ, ss_summ, '\t')){
//                        s++;
//                    }
//                    cout << tmp_ss << '\t' ;
//                    cout <<  mf1_dict[name][k]<<i << '\t' << s << '\t' << tmp_ss<< '\n';
                    if(mf1_dict[name][k].back()=='\t'){
                    mf1_dict[name][k]+=(tmp_ss);
                    }else{
                     mf1_dict[name][k]+='\t'+(tmp_ss);
                    }
                }
//                cout << mf1_dict[name].size() << endl;
            }
            ++i;
          }
        }    
       }else{
//       cout << mf1_dict["0"][0] << endl;
          stringstream ss(tmp);
          i = 0;
          while( getline(ss, tmp_ss, '\t') ){
            if( index.count(i) ){
               mf1_dict["0"][0] += '\t' + tmp_ss ;
//               cout << tmp_ss << '\t' << endl;
            }  
            ++i;          
          }        
       }     
       ++j;
    }
    cout << mf1_dict.size() << endl;
    ofstream mt_file ( mt_name );
    for( auto iter = mf1_dict.begin(); iter != mf1_dict.end(); ++iter ){
       for(string sub : iter->second ){
//            stringstream ss(sub);
//            while(getline(ss, tmp,'\t')){
//            cout << tmp << endl;
                mt_file << sub << '\n';
                
//            }
        }
    }
    mt_file.close();
//    
    cout << "Merge " << mf1_name << " & " << mf2_name << " to " << mt_name << "\n";
}
void remove_NA_col ( string txt_name, string to_name ){
   ifstream txt_file( txt_name );
    if(!txt_file) std::cerr << "Could not open " << txt_name << " file! \n";
    string tmp, tmp_ss;
    int i(0);
    ofstream to_file( to_name );
    while( getline(txt_file, tmp) ){
        if(tmp.substr(tmp.length()-2)=="NA"){
            continue;
        }else{
            to_file << tmp << endl;
        }      
    }
    cout << "Write " << txt_name << " to " << to_name << endl;
}
bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it) ) ++it;
    return !s.empty() && it == s.end();
}
void string_to_numeric( string txt_name, string ft_name, map<string, int> dict){
  ifstream txt_file( txt_name );
  if(!txt_file) std::cerr << "Could not open " << txt_name << " file! \n";
  int i(0), j(0);
  string tmp, tmp_ss, entity;
  ofstream ft_file( ft_name );
  while(getline(txt_file, tmp)){
    if( i > 0){
        stringstream ss(tmp);
        while(getline(ss, tmp_ss, '\t')){
            if(tmp_ss == "NA"){
               entity = "0";
            }else if(dict.count(tmp_ss)){
                entity = to_string(dict[tmp_ss]);
            }else if(isFloat(tmp_ss) || tmp_ss.substr(0,4)=="TCGA"){
//            cout << tmp_ss <<endl;
                entity = tmp_ss;
            }else {
                entity = "1";
            }
            ft_file << entity << '\t';
        }
    }else{
        ft_file << tmp;
    }
    ft_file << '\n';
    ++i;
  }
  txt_file.close();
  ft_file.close();
  cout << "Make " << txt_name << " numeric to " << ft_name << " done!" << endl; 
}
void remove_fc_fr(string txt_name, string ft_name, string fr_name, string fc_name){
    ifstream txt_file( txt_name );
  if(!txt_file) std::cerr << "Could not open " << txt_name << " file! \n";
  string tmp, tmp_ss;
  int i(0), j(0);
  read_first_row(txt_name, fr_name);
  read_first_col(txt_name, fc_name);
  vector<vector<string>> read_str;
  while(getline(txt_file, tmp)){
    stringstream ss(tmp);
     if(i > 0){
         j = 0;
         vector<string> tmp_vec;
         while(getline(ss, tmp_ss, '\t') ){
            if( j > 0 ){
                tmp_vec.emplace_back(tmp_ss);
            }
            ++j;
         }            
         read_str.emplace_back(tmp_vec);
     }   
    ++i;
  }
  txt_file.close();
  cout << read_str.size() << '\t' << read_str[0].size() << endl;

  ofstream ft_file( ft_name );
  for(int k = 0; k < read_str.size(); ++k){
    for(int e = 0; e < read_str[k].size(); ++e){
        if(e != read_str[k].size()-1){
            ft_file << read_str[k][e] << '\t';
        }else{
            ft_file << read_str[k][e];
        }
    }
    ft_file << endl;
  }
  ft_file.close();
  cout << "Remove first row and column from " << txt_name << " " << ft_name << " done! \n";
}
void count_num_row(string txt_name){
    ifstream txt_file( txt_name );
  if(!txt_file) std::cerr << "Could not open " << txt_name << " file! \n";
  string tmp, tmp_ss;
  int i(0), j(0);
  map<int, vector<int>> summary;
  while(getline(txt_file, tmp) ){
    stringstream ss(tmp);
    j = 0;
    while(getline(ss, tmp_ss, '\t')){
        ++j;
    }
    summary[j].emplace_back(i);
    ++i;
  }
  for(auto iter = summary.begin(); iter!=summary.end(); ++iter){
    cout << iter->first << " number of elements in  "<< iter->second.size() << " rows"<< endl;
//    if(iter->first==302){
//    for(int k = 0; k < iter->second.size(); ++k ){
//        cout <<  iter->second[k]<<" th row" << '\t';
//    }
//    cout << endl;
//    
//    }
  }
//  cout << "There are " << (j) << " elements in one row of " << txt_name << endl;
}
int main(){
  string tsv_name ("./DATA/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriverLoose_aggregated_matrix.tsv");
  string tsv_txt_name ( "data_mc3_gene.txt" );
  string xlsx_name ( "./DATA/1-s2.csv" );
  string xlsx_txt_name ( "data_1_s2_mmc.txt" );
  
//  readTSV_toTXT(xlsx_name, xlsx_txt_name);
  string symbol_name( "summary_mc3_gene.txt" );
  string xlsxsymbol_name( "summary_1_s2_mmc.txt" );
//  readTXT_symbols(xlsx_txt_name, xlsxsymbol_name);
  string fc_( "fc_mc3_gene.txt" );
  string fc_1( "fc_1_s2_mmc.txt" );
//  read_first_col(tsv_txt_name, fc_);
//  read_first_col(xlsx_txt_name, fc_1);
  string fr_( "fr_mc3_gene.txt" );
  string fr_1( "fr_1_s2_mmc.txt" );
//  read_first_row(tsv_txt_name, fr_);
//  read_first_row(xlsx_txt_name, fr_1);
  string tcga_code( "fr_mc3_gene_unique.txt" );
  string bar_name ( "fr_mc3_gene_unique_tcga.txt" );
  string ind_name ( "fr_mc3_gene_unique_ind.txt" ); 
  string summ_name ("fr_mc3_gene_unique_summary.txt");
//  read_TCGA( fr_, tcga_code, bar_name, ind_name, summ_name );
  string trans_symbol( "data_mc3_gene_transpose.txt" );
//  transpose_txt_file( tsv_txt_name, trans_symbol );
  string intersect_trans( "data_mc3_gene_transpose_1_s2_mmc.txt" );
//  intersection_first_col(fc_1, trans_symbol, intersect_trans);

  string fc_inter_trans( "fc_mc3_gene_transFc1_unique.txt" );
//  read_first_col(intersect_trans, fc_inter_trans);
  string tcga_code_( "fc_mc3_gene_transFc1_unique.txt" );
  string bar_name_ ( "fc_mc3_gene_transFc1_unique_tcga.txt" );
  string ind_name_ ( "fc_mc3_gene_transFc1_unique_ind.txt" ); 
  string summ_name_ ("fc_mc3_gene_transFc1_unique_summary.txt");
//  read_TCGA( fc_inter_trans, tcga_code_, bar_name_, ind_name_, summ_name_ );
  string mt_name ("data_mc3_gene_transpose_1_s2_mmc_TCGA_LF.txt");
  set<int> index{1, 4};
//  merge_two (intersect_trans, xlsx_txt_name, mt_name, index);
  set<int> idlf{4};
  string mtLF_name ("data_mc3_gene_transpose_1_s2_mmc_LF.txt");
//  merge_two (intersect_trans, xlsx_txt_name, mtLF_name, idlf);
  set<int> idtcga{1};
  string mttcga_name ("data_mc3_gene_transpose_1_s2_mmc_TCGA.txt");
//  merge_two (intersect_trans, xlsx_txt_name, mttcga_name, idtcga);
  string mt_nan_name ("data_mc3_gene_transpose_1_s2_mmc_TCGA_LF_rNAN.txt");
//  remove_NA_col( mt_name, mt_nan_name);
  map<string, int> dict_tcga = {{"ACC",1},{"BLCA",2},{"BRCA",3},{"CESC",4},{"CHOL",5},{"COAD",6},{"ESCA",7},{"GBM",8},{"HNSC",9},{"KICH",10},{"KIRC",11},{"KIRP",12},{"LGG",13},{"LIHC",14},{"LUAD",15},{"LUSC",16},{"MESO",17},{"OV",18},{"PAAD",19},{"PCPG",20},{"PRAD",21},{"READ",22},{"SARC",23},{"SKCM",24},{"STAD",25},{"TCGA Study",26},{"TGCT",27},{"THCA",28},{"UCEC",29},{"UCS",30}, {"UVM",31}};
  string fr_mt_nan_name ("fr_data_mc3_gene_transpose_1_s2_mmc_TCGA_LF_rNAN.txt");
//  read_first_row(mt_nan_name, fr_mt_nan_name);
  int col_index = 300;
  string dict_index_name ("TCGAStudy_data_mc3_gene_transpose_1_s2_mmc_TCGA_LF_rNAN.txt");
  read_col_index(mt_nan_name, dict_index_name, col_index);
  string mt_nan_num_name ("data_mc3_gene_transpose_1_s2_mmc_TCGA_LF_rNAN_num.txt");
// string_to_numeric( mt_nan_name, mt_nan_num_name, dict_tcga); // string to numeric
 
 string fea_label("data_gene_TCGA_LF.txt");
 string fr_fea_label_label("fr_data_gene_TCGA_LF.txt");
 string fc_fea_label_label("fc_data_gene_TCGA_LF.txt");
// remove_fc_fr(mt_nan_num_name, fea_label, fr_fea_label_label, fc_fea_label_label);
 string summ_fea_label("summary_data_gene_TCGA_LF.txt");
// readTXT_symbols(fea_label, summ_fea_label);
 string mt_nan_num_name_symb("summary_data_mc3_gene_transpose_1_s2_mmc_TCGA_LF_rNAN_num.txt");
// readTXT_symbols(mt_nan_num_name, mt_nan_num_name_symb);
//  return 1;
// count_num_row(fea_label);
// string fr_trans_symbol("fr_data_mc3_gene_transpose.txt");
// read_first_row(trans_symbol, fr_trans_symbol);
}
