#define _USE_MATH_DEFINES
#include<cmath>
#include <math.h>
#include<Rcpp.h>
#include <vector>
#include "aux_functions.h"
#include "tree.h"

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// The line above (depends) it will make all the depe ndcies be included on the file
using namespace Rcpp;


// void grow_fixed(Tree &new_tree,
//                 const Rcpp::NumericMatrix &x_train,
//                 const Rcpp::NumericMatrix &x_test,
//                 const Rcpp::NumericMatrix &xcut,
//
//                 int &n_min_size,
//                 int &id_t, // Object to identify if the same tree was returned
//                 int &verb_node_index // Store the node which was growned
// ){
//
//   int n_obs_train = x_train.rows();
//   int n_obs_test = x_test.rows();
//
//
//   int n_terminals = xcut.rows();
//   Rcpp::NumericVector new_node_train_index;
//   Rcpp::NumericVector new_node_test_index;
//
//   for(int i = 0; i<10;i++){
//   // Getting the new nodes
//   node new_node(1,new_node_train_index,
//                      new_node_test_index,
//                      -1,-1,1,
//                      0,xcut(i,0),0);
//
//   new_tree.list_node.push_back(new_node);
//   }
//
//   // Adding the node observations in the
//   for(int j = 0; j< 10; j++){
//
//     for(int i = 0 ; i< n_obs_train;i++){
//
//     if(x_train(i,0)<new_tree.list_node[j].var_split){
//
//     }
//
//     }
// }


void grow(Tree &new_tree,
          const Rcpp::NumericMatrix &x_train,
          const Rcpp::NumericMatrix &x_test,
          const Rcpp::NumericMatrix &xcut,

          int &n_min_size,
          int &id_t, // Object to identify if the same tree was returned
          int &verb_node_index, // Store the node which was growned
          Rcpp::NumericVector &change_grow_rules
){


  // Getting information from x
  int p = x_train.cols();
  int n_terminals;
  int n_nodes = new_tree.list_node.size();
  int g_node, g_var, n_train, n_test;
  double g_rule;


  // Getting observations that are on the left and the ones that are in the right
  Rcpp::NumericVector new_left_train_index;
  Rcpp::NumericVector new_right_train_index;
  Rcpp::NumericVector curr_obs_train; // Observations that belong to that terminal node

  // Same from above but for test observations
  // Getting observations that are on the left and the ones that are in the right
  Rcpp::NumericVector new_left_test_index;
  Rcpp::NumericVector new_right_test_index;
  Rcpp::NumericVector curr_obs_test; // Observations that belong to that terminal node


  // Creating the new splits for the xcut
  Rcpp::NumericVector new_left_xcut;
  Rcpp::NumericVector new_right_xcut;

  // Getting terminal nodes
  vector<node> t_nodes = new_tree.getTerminals();
  n_terminals = t_nodes.size();

  // Sampling one of the nodes (this function will sample the node index from
  // the vector)
  g_node = sample_int(n_terminals);
  g_var = sample_int(p);

  // Number of observations in the node that will grow
  curr_obs_train = t_nodes[g_node].obs_train; // Current observations on that node
  n_train = t_nodes[g_node].obs_train.size();

  // Doing the same but for the test
  curr_obs_test = t_nodes[g_node].obs_test; // Current observations on that node
  n_test = t_nodes[g_node].obs_test.size();

  // Do not let split if is not possible to split
  if( floor(n_train/2) <= n_min_size){
    id_t = 1;
    return;
  }

  // Selecting a split rule to choose
  // NEED TO SELECT ONLY SPLIT RULE FOR THAT TREE
  Rcpp::NumericVector x_curr(n_train);

  for(int i=0;i<n_train;i++){
    x_curr(i) = x_train(curr_obs_train(i),g_var);
  }
  // cout << " ERROR 1" << endl;
  // cout << " Size x_curr " << x_curr.size() << endl;

  // Sampling the rule
  g_rule = sample_rule(xcut(_,g_var),x_curr,n_min_size);
  // g_rule = sample_double(xcut(_,g_var));
  change_grow_rules.push_back(g_rule); // Adding this rule
  if(g_rule == 0.0003){
    id_t = 1;
    return;
  }

  // cout << "ERROR 2" << endl;

  // Iterating over the train observations
  for(int i = 0; i<n_train; i++){
    if(x_train(curr_obs_train[i],g_var)<=g_rule){
      new_left_train_index.push_back(curr_obs_train(i));
    } else {
      new_right_train_index.push_back(curr_obs_train(i));
    }
  }

  // Iterating over the test observations
  for(int i = 0; i<n_test; i++){
    if(x_test(curr_obs_test[i],g_var)<=g_rule){
      new_left_test_index.push_back(curr_obs_test[i]);
    } else {
      new_right_test_index.push_back(curr_obs_test[i]);
    }
  }

  // cout << "Right node size " << new_right_train_index.size() << endl;
  // cout << "Left node size " << new_left_train_index.size() << endl;
  // cout << "That boolean is true " << ((new_right_train_index.size()>=n_min_size) & (new_left_train_index.size()>=n_min_size)) << endl;

  if(((new_right_train_index.size()>=n_min_size) & (new_left_train_index.size()>=n_min_size))){

    // cout << "YESSS" << endl;
    // Getting the new nodes
    node new_node_left(n_nodes,new_left_train_index,
                       new_left_test_index,
                       -1,-1,t_nodes[g_node].depth+1,
                       g_var,g_rule,t_nodes[g_node].mu);

    node new_node_right(n_nodes+1,new_right_train_index,
                        new_right_test_index,
                        -1,-1,t_nodes[g_node].depth+1,
                        g_var,g_rule,t_nodes[g_node].mu);

    // Adding the new nodes
    new_tree.list_node.push_back(new_node_left);
    new_tree.list_node.push_back(new_node_right);

    // Transforming the growned node into non-terminal
    new_tree.list_node[t_nodes[g_node].index].left = n_nodes;
    new_tree.list_node[t_nodes[g_node].index].right = n_nodes+1;

    // Storing the index of the node index
    verb_node_index = t_nodes[g_node].index;
  }

  return;


}




// Prune a tree
void prune(Tree &new_tree, int& id_t, int& verb_node_index){

  // New list of nodes
  int n_nodes = new_tree.list_node.size();
  vector<node> new_nodes;

  // Can't prune a root
  if(new_tree.list_node.size()==1){
    id_t = 1;
    return;
  }

  // Getting the parents of terminal nodes
  vector<node> parent_left_right;
  for(int i=0;i<n_nodes;i++){
    if(!new_tree.list_node[i].isTerminal()){
      if(new_tree.list_node[new_tree.list_node[i].left].isTerminal()){ // IS THE INDEX NOT THE POSITION!
        if(new_tree.list_node[new_tree.list_node[i].right].isTerminal()){
          parent_left_right.push_back(new_tree.list_node[i]); // Adding the parent
        }
      }
    }
  }


  // Getting the number of internal and the node to be pruned
  int n_parent = parent_left_right.size();
  int p_node = sample_int(n_parent);


  // Returning to be a terminal node
  int left_index = new_tree.list_node[parent_left_right[p_node].index].left;
  int right_index = new_tree.list_node[parent_left_right[p_node].index].right;
  new_tree.list_node[parent_left_right[p_node].index].left=-1;
  new_tree.list_node[parent_left_right[p_node].index].right=-1;

  // Storing the index of the node that was pruned
  verb_node_index = parent_left_right[p_node].index;

  // Adding new trees

  for(int k = 0;k<new_tree.list_node.size();k++){
    if((new_tree.list_node[k].index != left_index) && (new_tree.list_node[k].index != right_index)){


      // Doing for trees greater than left index
      if(new_tree.list_node[k].index>left_index){ /// COULD ALSO USE IF new_tree.list_node[k].depth >depth_prune_node
        new_tree.list_node[k].index = new_tree.list_node[k].index-2;
      }
      // Checking if the left index is greater than it should be
      if(new_tree.list_node[k].left>left_index){
        new_tree.list_node[k].left = new_tree.list_node[k].left-2;
      }
      // Checking if the right index is greater than it should be
      if(new_tree.list_node[k].right>left_index){
        new_tree.list_node[k].right = new_tree.list_node[k].right-2;
      }

      new_nodes.push_back(new_tree.list_node[k]);
    }
  }

  // Updating the new nodes
  new_tree.list_node = new_nodes;
  return ;
}

// Change a tree
void change(Tree &new_tree,
            const Rcpp::NumericMatrix x_train,
            const Rcpp::NumericMatrix x_test,
            const Rcpp::NumericMatrix xcut,
            int n_min_size,
            int &id_t,
            int &verb_node_index,
            Rcpp::NumericVector change_grow_rules){

  // Declaring important values and variables
  int n_nodes = new_tree.list_node.size();
  int p = x_train.cols();
  int n_train,n_test;
  int c_var,c_node;
  double c_rule;

    // Getting observations that are on the left and the ones that are in the right
  Rcpp::NumericVector new_left_train_index;
  Rcpp::NumericVector new_right_train_index;
  Rcpp::NumericVector curr_obs_train; // Observations that belong to that terminal node

  // Getting observations that are on the left and the ones that are in the right
  Rcpp::NumericVector new_left_test_index;
  Rcpp::NumericVector new_right_test_index;
  Rcpp::NumericVector curr_obs_test; // Observations that belong to that terminal node


  // Can't prune a root
  if(new_tree.list_node.size()==1){
    id_t = 1;
    return;
  }

  // Getting the parents of terminal nodes (NOG)
  vector<node> parent_left_right;
  for(int i=0;i<n_nodes;i++){
    if(new_tree.list_node[i].isTerminal()==0){
      if(new_tree.list_node[new_tree.list_node[i].left].isTerminal()==1 && new_tree.list_node[new_tree.list_node[i].right].isTerminal()==1){
        parent_left_right.push_back(new_tree.list_node[i]); // Adding the parent
      }
    }
  }


  // Getting the number of parent of terminals and selecting the one to be changed
  int n_parent = parent_left_right.size();

  // Selecting the node and var
  c_node = sample_int(n_parent);
  c_var = sample_int(p);



  // Number of train observations that will be changed
  curr_obs_train = parent_left_right[c_node].obs_train;
  n_train = parent_left_right[c_node].obs_train.size();

  // Number of test observations that will be changed
  curr_obs_test = parent_left_right[c_node].obs_test;
  n_test = parent_left_right[c_node].obs_test.size();

  // NEED TO SELECT ONLY SPLIT RULE FOR THAT TREE
  Rcpp::NumericVector x_curr(n_train);

  for(int i=0;i<n_train;i++){
    x_curr(i) = x_train(curr_obs_train(i),c_var);
  }

  // Getting the c_rule
  c_rule = sample_rule(xcut(_,c_var),x_curr,n_min_size);
  // c_rule = sample_double(xcut(_,c_var));


  change_grow_rules.push_back(c_rule); // Adding this rule

  if(c_rule == 0.0003){
       id_t = 1;
    return;
  }

  // Iterating over the train observations
  for(int i = 0; i<n_train; i++){
    if(x_train(curr_obs_train[i],c_var)<=c_rule){
      new_left_train_index.push_back(curr_obs_train[i]);
    } else {
      new_right_train_index.push_back(curr_obs_train[i]);
    }
  }

  // Iterating over the test observations
  for(int i = 0; i<n_test; i++){
    if(x_test(curr_obs_test[i],c_var)<=c_rule){
      new_left_test_index.push_back(curr_obs_test[i]);
    } else {
      new_right_test_index.push_back(curr_obs_test[i]);
    }
  }


  // Verifying if is the correct min node size
  if((new_right_train_index.size()>=n_min_size) & (new_left_train_index.size()>=n_min_size)){

    // Returning the nodes that will be changed
    int left_index = new_tree.list_node[parent_left_right[c_node].index].left;
    int right_index = new_tree.list_node[parent_left_right[c_node].index].right;


    // Changing the rules and observations that belong to each one of
    // that nodes
    new_tree.list_node[left_index].var = c_var;
    new_tree.list_node[right_index].var = c_var;
    new_tree.list_node[left_index].var_split = c_rule;
    new_tree.list_node[right_index].var_split = c_rule;
    new_tree.list_node[left_index].obs_train = new_left_train_index;
    new_tree.list_node[right_index].obs_train = new_right_train_index;


    // Changing the test observations
    new_tree.list_node[left_index].obs_test = new_left_test_index;
    new_tree.list_node[right_index].obs_test = new_right_test_index;

    // storing the index of the changed node
    verb_node_index = parent_left_right[c_node].index;

  } else {
    id_t = 1;
    verb_node_index = parent_left_right[c_node].index;
  }
  return;

};

double node_loglikelihood(Rcpp::NumericVector residuals,
                          node current_node,
                          double tau,
                          double tau_mu) {

  // Decarling the quantities
  int n_size = current_node.obs_train.size();
  Rcpp::NumericVector current_node_index = current_node.obs_train;
  double sum_sq_r = 0 , sum_r = 0;

  for(int i = 0;i<n_size;i++){
    sum_sq_r+=residuals(current_node_index(i))*residuals(current_node_index(i));
    sum_r+=residuals(current_node_index(i));

  }

  return -0.5*tau*sum_sq_r+0.5*((tau*tau)*(sum_r*sum_r))/(tau*n_size+tau_mu)-0.5*log(tau*n_size+tau_mu);
}

//[[Rcpp::export]]
double node_loglikelihood_r(Rcpp::NumericVector residuals,
                            Rcpp::NumericVector node_index,
                          double tau,
                          double tau_mu) {

  // Decarling the quantities
  int n_size = residuals.size();
  double sum_sq_r = 0 , sum_r = 0;

  for(int i = 0;i<n_size;i++){
    sum_sq_r+=residuals(i)*residuals(i);
    sum_r+=residuals(i);

  }

  return -0.5*tau*sum_sq_r+0.5*((tau*tau)*(sum_r*sum_r))/(tau*n_size+tau_mu)-0.5*log(tau*n_size+tau_mu);
}

double tree_loglikelihood(Rcpp::NumericVector residuals,
                          Tree curr_tree,
                          double tau,
                          double tau_mu
) {

  // Declaring important quantities
  vector<node> terminal_nodes = curr_tree.getTerminals();
  int number_nodes = terminal_nodes.size();
  double loglike_sum = 0;


  for(int i = 0; i<number_nodes; i++) {
    loglike_sum+=node_loglikelihood(residuals,
                                    terminal_nodes[i],tau,tau_mu);
  }

  return loglike_sum;
}

double tree_loglikelihood_verb(Rcpp::NumericVector residuals,
                               Tree new_tree,
                               Tree curr_tree,
                               double verb,
                               int verb_node_index,
                               double tau,
                               double tau_mu){
  if(verb < 0.28){ // Grow prob
    return (node_loglikelihood(residuals,new_tree.list_node[new_tree.list_node[verb_node_index].left],tau,tau_mu)+
            node_loglikelihood(residuals,new_tree.list_node[new_tree.list_node[verb_node_index].right],tau,tau_mu)) - node_loglikelihood(residuals,new_tree.list_node[verb_node_index],tau,tau_mu);
  } else if( verb <= 0.56){ // Prune prob
    return node_loglikelihood(residuals, curr_tree.list_node[verb_node_index],tau,tau_mu)-(node_loglikelihood(residuals, curr_tree.list_node[curr_tree.list_node[verb_node_index].left],tau,tau_mu)+node_loglikelihood(residuals, curr_tree.list_node[curr_tree.list_node[verb_node_index].right],tau,tau_mu));
  } else { // Change case
    return (node_loglikelihood(residuals,new_tree.list_node[new_tree.list_node[verb_node_index].left],tau,tau_mu)+node_loglikelihood(residuals,new_tree.list_node[new_tree.list_node[verb_node_index].right],tau,tau_mu))-(node_loglikelihood(residuals, curr_tree.list_node[curr_tree.list_node[verb_node_index].left],tau,tau_mu)+node_loglikelihood(residuals, curr_tree.list_node[curr_tree.list_node[verb_node_index].right],tau,tau_mu));
  }

}



// Updating the \mu
Tree update_mu(Rcpp::NumericVector residuals,
               Tree curr_tree,
               double tau,
               double tau_mu){

  // Declaring important quantities
  vector<node> terminal_nodes = curr_tree.getTerminals();
  int number_nodes = terminal_nodes.size();
  double mu_mean, mu_sd;
  double sum_residuals;
  int nj = 0;
  // Iterating over terminal nodes
  for(int i = 0;i<number_nodes;i++){
    sum_residuals = 0;
    nj = terminal_nodes[i].obs_train.size();

    // Calculating the sum of residuals for that terminal node
    for(int j = 0;j<nj;j++){
      sum_residuals+=residuals(terminal_nodes[i].obs_train[j]);
    }
    mu_mean = (tau*sum_residuals)/(nj*tau+tau_mu);
    // cout << "Number obs node terminal"<< nj << endl;
    mu_sd = sqrt(1/(nj*tau+tau_mu));
    curr_tree.list_node[terminal_nodes[i].index].mu = R::rnorm( mu_mean, mu_sd);
    // curr_tree.list_node[terminal_nodes[i].index].mu = mu_mean;


  }
  // cout << "Mu mean: " << mu_mean << endl;
  // cout << "Mu sd: " << mu_sd<< endl;
  return curr_tree;
}

// Calculate the density of a half cauchy location 0 and sigma
//[[Rcpp::export]]
double dhcauchy(double x, double location, double sigma){

  if( x>location) {
    return (1/(M_PI_2*sigma))*(1/(1+((x-location)*(x-location))/(sigma*sigma)));
  } else {
    return 0.0;
  }
}


double update_tau_old(Rcpp::NumericVector y,
                      Rcpp::NumericVector y_hat,
                      double a_tau,
                      double d_tau){

  // Function used in the development of the package where I checked
  // contain_nan(y_hat);
  int n = y.size();
  double sum_sq_res = 0;
  for(int i = 0;i<n;i++){
    sum_sq_res = sum_sq_res + (y(i)-y_hat(i))*(y(i)-y_hat(i));
  }
  return R::rgamma((0.5*n+a_tau),1/(0.5*sum_sq_res+d_tau));
}

// double update_tau(arma::vec y,
//                   arma::vec y_hat,
//                   double naive_sigma,
//                   double curr_tau){
//
//   int n=y.size();
//   double curr_sigma, proposal_tau,proposal_sigma, acceptance;
//
//   curr_sigma = 1/(sqrt(curr_tau));
//   // Calculating the residuals
//   arma::vec res = (y-y_hat);
//
//   // Getting the sum squared of residuals
//   double sum_sq_res=dot(res,res);
//
//   // Calculating manually
//   // double sum_sq_res = 0;
//   // for(int i=0;i<y.size();i++){
//   //   sum_sq_res = sum_sq_res + (y(i)-y_hat(i))*(y(i)-y_hat(i));
//   // }
//
//
//   // Getting a proposal sigma
//   proposal_tau = R::rgamma(n/2+1,2/(sum_sq_res));
//   // proposal_tau = R::rgamma(51,1.76);
//
//   proposal_sigma  = 1/(sqrt(proposal_tau));
//
//   acceptance = exp(log(dhcauchy(proposal_sigma,0,naive_sigma))+3*log(proposal_sigma)-log(dhcauchy(curr_sigma,0,naive_sigma))-3*log(curr_sigma));
//   // return proposal_tau;
//   // return proposal_tau;
//   if(R::runif(0,1)<=acceptance){
//     return proposal_tau;
//   } else{
//     return curr_tau;
//   }
// }
//
//
//
//
//
void get_prediction_tree(Tree curr_tree,
                         const Rcpp::NumericMatrix x_train,
                         const Rcpp::NumericMatrix x_test,
                         Rcpp::NumericVector &prediction_train,
                         Rcpp::NumericVector &prediction_test){



  // Getting terminal nodes
  vector<node> terminal_nodes = curr_tree.getTerminals();
  int n_terminals = terminal_nodes.size();


  vector<int> curr_obs_train;
  vector<int> new_left_train_index, new_right_test_index;

  // Iterating to get the predictions for a training test
  for(int i = 0;i<n_terminals;i++){

    // Updating the training predictions
    for(int j = 0; j<terminal_nodes[i].obs_train.size();j++){
      prediction_train[terminal_nodes[i].obs_train[j]] = terminal_nodes[i].mu;
    }

    for(int t = 0; t<terminal_nodes[i].obs_test.size();t++){
      prediction_test[terminal_nodes[i].obs_test[t]] = terminal_nodes[i].mu;
    }

  }

}


double tree_log_prior_verb(Tree new_tree,
                           Tree curr_tree,
                           double verb,
                           double alpha,
                           double beta) {

  double log_tree_p = 0.0;

  if (verb<0.56){

    // Adding for new tree
    for(int i=0;i<curr_tree.list_node.size();i++){
      if(curr_tree.list_node[i].isTerminal()){
        log_tree_p = log_tree_p - log(1-alpha/pow((1+curr_tree.list_node[i].depth),beta));
      } else {
        log_tree_p = log_tree_p - (log(alpha)-beta*log(1+curr_tree.list_node[i].depth));
      }
    }

    // Subtracting for the current tree
    for(int i=0;i<new_tree.list_node.size();i++){
      if(new_tree.list_node[i].isTerminal()){
        log_tree_p = log_tree_p + log(1-alpha/pow((1+new_tree.list_node[i].depth),beta));
      } else {
        log_tree_p = log_tree_p + (log(alpha)-beta*log(1+new_tree.list_node[i].depth));
      }
    }

  }

  return log_tree_p;
}

// Get the log transition probability
double log_transition_prob(Tree curr_tree,
                           Tree new_tree,
                           double verb){

  // Getting the probability
  double log_prob = 0;
  // In case of Grow: (Prob from Grew to Current)/(Current to Grow)
  if(verb < 0.3){
    log_prob = log(0.3/new_tree.n_nog())-log(0.3/curr_tree.n_terminal());
  } else if (verb < 0.6) { // In case of Prune: (Prob from the Pruned to the current)/(Prob to the current to the prune)
    log_prob = log(0.3/new_tree.n_terminal())-log(0.3/curr_tree.n_nog());
  }; // In case of change log_prob = 0; it's already the actual value

  return log_prob;

}

//[[Rcpp::export]]
List bart(const Rcpp::NumericMatrix x_train,
          const Rcpp::NumericVector y,
          const Rcpp::NumericMatrix x_test,
          const Rcpp::NumericMatrix xcut,
          int n_tree,
          int n_mcmc,
          int n_burn,
          int n_min_size,
          double tau, double mu,
          double tau_mu, double naive_sigma,
          double alpha, double beta,
          double a_tau, double d_tau){

  // Declaring common variables
  double verb;
  double acceptance;
  double log_transition_prob_obj;
  double acceptance_ratio = 0;
  int post_counter = 0;
  int id_t ,verb_node_index; // Id_t: Boolean to verify if the same tree was generated
  // Verb_node_index: Explicit the node the was used;
  Rcpp::NumericVector n_nodes;
  // Getting the number of observations
  int n_train = x_train.rows();
  int n_test = x_test.rows();
  Rcpp:: NumericVector verb_list;
  int id_t_c = 0;

  // Creating a vector to store all split rules that are being proposed
  Rcpp::NumericVector change_grow_rules;

  // Creating the variables
  int n_post = (n_mcmc-n_burn);
  Rcpp::NumericMatrix y_train_hat_post(n_post,n_train);
  Rcpp::NumericMatrix y_test_hat_post(n_post,n_train);

  Rcpp::NumericVector tau_post;

  Tree new_tree(n_train,n_test);

  // Creating a vector of multiple trees
  vector<Tree> current_trees;
  for(int i = 0; i<n_tree;i++){
    current_trees.push_back(new_tree);

  }


  // Creating a matrix of zeros of y_hat
  y_train_hat_post.fill(0);
  y_test_hat_post.fill(0);

  // Creating the partial residuals and partial predictions
  Rcpp::NumericVector partial_pred(n_train),partial_residuals(n_train);
  Rcpp::NumericVector prediction_train(n_train),prediction_test(n_test); // Creating the vector only one prediction
  Rcpp::NumericVector prediction_test_sum(n_test);

  // Initializing zero values
  partial_pred.fill(0);
  partial_residuals.fill(0);
  prediction_train.fill(0);
  prediction_test.fill(0);


  // Iterating over all MCMC samples
  for(int i=0;i<n_mcmc;i++){

    prediction_test_sum.fill(0);


    // Iterating over the trees
    for(int t = 0; t<n_tree;t++){

      // Initializing if store tree or not
      id_t = 0;
      // Rcout << std::hex << "END "<<i<<endl;
      // Updating the prediction tree and prediction test
      // cout << "Pred(Before) : " << prediction_train(0) << endl;
      // get_prediction_tree(current_trees[t],x_train,x_test,prediction_train,prediction_test);
      // cout << "Pred(After) : " << prediction_train(0) << endl;

      // Updating the residuals
      partial_residuals = y ;// - (partial_pred - prediction_train);

      // Sampling a verb. In this case I will sample a double between 0 and 1
      // Grow: 0-0.3
      // Prune 0.3 - 0.6
      // Change: 0.6 - 1.0
      // Swap: Not in this current implementation
      verb = R::runif(0,1);
      // Forcing the first trees to grow
      if(current_trees[t].list_node.size()==1 ){
        verb = 0.1;
      }

      // Proposing a new tree
      // if(verb<=0.5){
      //   grow(new_tree, x_train,x_test,xcut, n_min_size,id_t,verb_node_index,change_grow_rules);
      // } else if ( verb>0.4 & verb <= 0.7){
      //   prune(new_tree,id_t,verb_node_index);
      //   // }
      // } else if (verb>0.7 & verb <= 1){
      //   change(new_tree,x_train,x_test,xcut,n_min_size,id_t,verb_node_index,change_grow_rules);
      // }

      if(verb<=0.5){
        grow(new_tree, x_train,x_test,xcut, n_min_size,id_t,verb_node_index,change_grow_rules);
      } else {
        prune(new_tree,id_t,verb_node_index);
      }


      // Calculating or not the likelihood (1 is for case where the trees are the same)
      // if(id_t == 0){

        // No new tree is proposed (Jump log calculations)
        if( (verb <=0.7) && (current_trees[t].list_node.size()==new_tree.list_node.size())){
          log_transition_prob_obj = 0;
        } else {
          log_transition_prob_obj = log_transition_prob(current_trees[t],new_tree,verb);

        }

        // acceptance = tree_loglikelihood_verb(partial_residuals,new_tree,current_trees[t],verb,verb_node_index,tau,tau_mu) + tree_log_prior_verb(new_tree,current_trees[t],verb,alpha,beta)+ log_transition_prob_obj;
        acceptance = tree_loglikelihood(partial_residuals,new_tree,tau,tau_mu) - tree_loglikelihood(partial_residuals,current_trees[t],tau,tau_mu) + tree_log_prior_verb(new_tree,current_trees[t],verb,alpha,beta) + log_transition_prob_obj;

        // cout << "Prob. Acceptance" << exp(acceptance) << endl;

        // Testing if will acceptance or not
        if( (R::runif(0,1)) < exp(acceptance)){
          acceptance_ratio++;
          // cout << "ACCEPTED" << endl;
          current_trees[t] = new_tree;
        }

      // } // End the test likelihood calculating for same trees
      if(id_t==1){
        id_t_c++;
      }

      // Generating new \mu values for the accepted (or not tree)
      current_trees[t] = update_mu(partial_residuals,current_trees[t],tau,tau_mu);

      // Updating the prediction train and prediction test with the sampled \mu values
      get_prediction_tree(current_trees[t],x_train,x_test,prediction_train,prediction_test);

      // partial_pred = y + prediction_train - partial_residuals;
      partial_pred = prediction_train ;

      // cout << "Tree size" << t <<" : " << current_trees[t].list_node.size() << endl;

      prediction_test_sum += prediction_test;
    }

    // Updating tau
    // past_tau = tau;
    // cout << "Y(0)" << (0) << endl;
    tau = update_tau_old(y,partial_pred,a_tau,d_tau);
    // tau = update_tau_old(y,partial_pred,a_tau,d_tau);
    // cout << "TAU: " << tau << endl;
    n_nodes.push_back(current_trees[0].list_node.size());

    // Updating the posterior matrix
    if(i>=n_burn){
      y_train_hat_post(post_counter,_) = partial_pred;
      y_test_hat_post(post_counter,_) = prediction_test_sum;
      //Updating tau
      tau_post.push_back(tau);
      post_counter++;
    }

  }
  cout << "Acceptance Ratio = " << acceptance_ratio/n_tree << endl;
  cout << "Identitical Tree Counter= " << id_t_c << endl;

  return Rcpp::List::create(_["y_train_hat_post"] = y_train_hat_post,
                            _["y_test_hat_post"] = y_test_hat_post,
                            _["tau_post"] = tau_post,
                            _["n_nodes"] = n_nodes,
                            _["change_grow_rules"] = change_grow_rules);

}

//
// //[[Rcpp::export]]
// arma::vec sub_setting(arma::vec vector, int beg, int end){
//
//   return vector.subvec(beg,end);
// }

// Small tests to find the bug on my implementation of cpp bart
// //[[Rcpp::export]]
// double bart_tests(arma::mat x,
//                   arma::vec y,
//                   arma::mat x_test,
//                   double tau,
//                   double tau_mu){
//
//   // Initializing a tree node
//   Tree simple_tree(x.n_rows,x_test.n_rows);
//   // return tree_loglikelihood(y,simple_tree,tau,tau_mu);
//   simple_tree = update_mu(y,simple_tree,tau,tau_mu);
//   return 0;
// }

// //[[Rcpp::export]]
// int initialize_test(Eigen::MatrixXd x_train,Eigen::MatrixXd x_test,
//                     Eigen::VectorXd residuals,
//                     double tau, double tau_mu, int node_min_size){
//
//   int id_t = 0;
//   int verb_node_index = -1;
//   int n_obs_train(x_train.rows()), n_obs_test(x_test.rows());
//   Tree tree1(n_obs_train,n_obs_test);
//   Tree new_tree = grow(tree1,x_train,x_test,node_min_size, &id_t,&verb_node_index);
//   id_t = 0;
//   verb_node_index = -1;
//   Tree new_tree_two = grow(new_tree,x_train,x_test,node_min_size, &id_t,&verb_node_index);
//   Tree new_tree_three = grow(new_tree_two,x_train,x_test,node_min_size, &id_t,&verb_node_index);
//   Tree new_tree_four = grow(new_tree_three,x_train,x_test,node_min_size, &id_t,&verb_node_index);
//   Tree new_tree_five = grow(new_tree_four,x_train,x_test,node_min_size, &id_t,&verb_node_index);
//   Tree new_tree_six = grow(new_tree_five,x_train,x_test,node_min_size, &id_t,&verb_node_index);
//
//
//   id_t = 0;
//   verb_node_index = -1;
//
//   new_tree_six.DisplayNodes();
//   cout << "======= SWAPP LINE =======" << endl;
//
//   Tree swap_tree = swap_old(new_tree_six,x_train,node_min_size);
//   // Tree swap_tree = prune(new_tree_five, &id_t,&verb_node_index);
//
//   swap_tree.DisplayNodes();
//   cout << " ==== " << endl;
//   cout << " Same Tree was generated: "<<id_t << endl;
//   cout << " ==== " << endl;
//
//   cout << " ==== " << endl;
//   cout << " The modified node was: "<< verb_node_index << endl;
//   cout << " ==== " << endl;
//
//
//   // Tree new_tree_three = grow(new_tree_two,x,2);
//   // Tree new_tree_four = grow(new_tree_three,x,2);
//   // Tree new_tree_five = grow(new_tree_four,x,2);
//   // Tree new_tree_six = grow(new_tree_five,x,2);
//   // new_tree_four.DisplayNodes();
//   // // Tree new_three_tree = grow(new_tree_two,x,2);
//
//   // cout << "======= SWAPP LINE =======" << endl;
//   // Tree swap_tree(n_obs);
//   // get_prediction_tree(new_tree_five,x,true);
//   // for(int k = 0;k<1000;k++){
//     // cout << k << endl;
//     // swap_tree = swap(new_tree_four,x,2);
//     // swap_tree.DisplayNodes();
//   // }
//   // Tree update_mu_tree = update_mu(residuals,change_tree_two,tau,tau_mu);
//   // Eigen::VectorXd pred_vec  = get_prediction_tree(update_mu_tree,x_new,false);
//   // // get_prediction_tree(change_tree_two,x,true);
//   // cout << "Pred x: ";
//   // for(int i = 0; i<pred_vec.size();i++){
//   //   cout  << pred_vec.coeff(i)<< " " ;
//   // }
//
//   // cout << "Tree prior" << tree_log_prior(new_tree_two,0.95,2) << endl;
//   // vector<int> p_index = new_tree_two.getParentTerminals();
//
//   // Tree prune_tree = prune(new_three_tree);
//   // prune_tree.DisplayNodes();
//
//   return 0;
//
// }

