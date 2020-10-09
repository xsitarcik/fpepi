#!/usr/bin/env python
# coding: utf-8

import random
from scipy import stats
import numpy as np
import pandas as pd
import os
from math import exp, log, gamma, sin,pi
from math import pow as mpow
from itertools import combinations #to generate combinations
from math import log2
from scipy.stats import entropy
from scipy.stats import chisquare
from scipy.stats import chi2

#######################################################################
#PLEASE SET THE INPUT PARAMETERS HERE
gini_population = 25
k2_population = 25
iteration_size = 50
alpha_value = 0.1
path = "/tf/70.1600.0.antesnp100.txt"
output_file = "results.txt"

#PLEASE SET THE HYPERPARAMETERS HERE (OPTIONAL)
init_prob = 0.2
min_dof = 5
index_beta = 1.5
required_iterations_for_taboo = 5
zeta_radius = 1.001
bee_limit = 5
reflection_coef = 1
extension_coef = 2
compression_coef = 0.5
shrink_coef = 0.5
elite_crossover = 0.5
######################################################################

#inprogram parameters
dimension = 2
lowest_coord = 0 #set later
highest_coord = 0 #set later
total_combinations = 0 #set later
number_of_SNPs = 0 #set later
best_solutions = []
population_size = gini_population + k2_population
num_of_samples = 0 #be set later

def compute_mi(pos):
    a = calc_joint_entropy(pos)
    b = calc_joint_entropy(-1)
    r = []
    for i in pos:
        r.append(i)
    r.append(-1)
    c = calc_joint_entropy(r)
    return (1-(a + b - c))

def calc_joint_entropy(pos):
    a,cts = np.unique(df.iloc[:,pos],return_counts=True,axis=0)
    cts = cts/num_of_samples
    return -sum(cts*np.log2(cts))

def get_scores(snp):
    global visited
    if(gini_computed[snp[0]][snp[1]])==-1:
        contingency_table = pd.crosstab(df.Class, [df.iloc[:,snp[0]],df.iloc[:,snp[1]]])
        gini_score = compute_gini(contingency_table)
        gini_computed[snp[0]][snp[1]] = gini_score
        k2_score = compute_K2_score(contingency_table)
        k2score_computed[snp[0]][snp[1]] = k2_score
        return gini_score,k2_score
    else:
        return gini_computed[snp[0]][snp[1]],k2score_computed[snp[0]][snp[1]]

def compute_gini(cont_table):
    Pi, gini = 0,0
    for j, k in cont_table:
        pi0 = cont_table[j][k][0] / num_of_samples;
        pi1 = cont_table[j][k][1] / num_of_samples;
        Pi = pi0 + pi1
        gini += Pi * (1 - ((pi0*pi0) + (pi1*pi1)))
    return gini

def compute_K2_score(cont_table):
    k2_score = 0
    for j, k in cont_table:
        sum_controls, sum_cases, sum_all = 0, 0, 0
        number_of_controls, number_of_cases = 0, 0
        for number_of_controls in np.arange(1,cont_table[j][k][0]+1):
            sum_controls= sum_controls + log(number_of_controls)
        for number_of_cases in np.arange(1,cont_table[j][k][1]+1):
            sum_cases = sum_cases + log(number_of_cases)
        for i in np.arange(number_of_cases+1,number_of_cases+number_of_controls+2):
            sum_all = sum_all + log(i)
        k2_score = k2_score + sum_all - sum_controls
    return k2_score
   
def is_the_same(coordinate):
    return True if coordinate[0]==coordinate[1] else False
    
def get_p_value(position):
    if precomputed_p_values[position[0]][position[1]] == -1:
        contingency_table = pd.crosstab(df.Class, [df.iloc[:,position[0]], df.iloc[:,position[1]]])
        for j, k in contingency_table:
            suma = contingency_table[j][k][0]+contingency_table[j][k][1]
            if suma<min_dof:
                contingency_table.drop((j, k), axis = 1,inplace=True)
        chi2_stat, p_val, dof, ex = stats.chi2_contingency(contingency_table,lambda_="log-likelihood")
        precomputed_p_values[position[0]][position[1]] = p_val
        return p_val
    else:
        return precomputed_p_values[position[0]][position[1]]   

def is_in_taboo_coord(coord,taboo):
    if coord[0]==taboo[0] or coord[0]==taboo[1] or coord[1]==taboo[0] or coord[1]==taboo[1]:
        return True
    else:
        return False
        
def generate_population(population_size,dimension):
    pop = np.random.randint(lowest_coord, highest_coord,(gini_population,dimension))
    for i,agent in enumerate(pop):
        while is_the_same(agent):
            pop[i] = np.random.randint(lowest_coord, highest_coord,dimension)
        pop[i] = standardize_coord(agent) 
    return pop

def round_coord(coordinate):
    snp = []
    for i in range(dimension):
        if coordinate[i]<=lowest_coord:
            snp.append(int(round(lowest_coord)))
        elif coordinate[i]>=highest_coord:
            snp.append(int(round(highest_coord)))
        else:
            snp.append(int(round(coordinate[i])))
    return snp

# dynamic switching probability strategy
def get_prob_switch(init_prob, num_of_iters, curr_iter):
    return init_prob - (0.1 * ((num_of_iters - curr_iter) / num_of_iters))

def standardize_coord(position):
    if position[0]>position[1]:
        aux = position[0]
        position[0] = position[1]
        position[1] = aux
    return position

def handle_same_SNPs(position):
    position = round_coord(position)
    if is_the_same(position):
        random_dim = random.randint(0,dimension-1)
        if position[random_dim]>=98.5:
            random_dir = -1
        elif position[random_dim]<0.5:
            random_dir = 1
        else: 
            random_dir = random.choice([-1,1])
        position[random_dim] += random_dir
    return standardize_coord(round_coord(position))


def levy_flight(index_beta):
    exponent = 1 / index_beta
    numerator = gamma(1 + index_beta) * sin((pi * index_beta) / 2)
    denominator = gamma((1 + index_beta) / 2) * index_beta * mpow(2, (index_beta - 1) / 2)
    base = mpow(numerator / denominator, exponent)
    base = mpow(base, 2)

    u = np.random.normal(0, base)
    v = abs(np.random.normal(0, 1))

    random_step = u / mpow(v, 1 / index_beta)

    return random_step

def local_search_elite(prev_flower, prev_random_flower_one, prev_random_flower_two,best):
    uni_dist = np.random.uniform(0, 1)
    one_two = [uni_dist*(i-j) for i,j in zip(prev_random_flower_one,prev_random_flower_two)]
    
    uni_dist = np.random.uniform(0, 1)
    best_dist = [uni_dist*(i-j) for i,j in zip(best,prev_flower)]
    
    suma = [i+j for i,j in zip(best_dist,one_two)]
    new_flower_tmp = [i+j for i,j in zip(prev_flower,suma)]
    return handle_same_SNPs(clip_agent(new_flower_tmp))

#def local_search(prev_flower, prev_random_flower_one, prev_random_flower_two):
#    prev = np.array(prev_flower)
#    prev_random_one = np.array(prev_random_flower_one)
#    prev_random_two = np.array(prev_random_flower_two)
#    uni_dist = np.random.uniform(0, 1)
#    prev_random_merged = prev_random_flower_one-prev_random_flower_two
#    random_local = [uni_dist * i for i in prev_random_merged]
#    new_flower_tmp = np.add(prev, random_local)
#    return handle_same_SNPs(clip_agent(new_flower_tmp))

def better_gini(new,old):
    return True if get_scores(new)[0]<get_scores(old)[0] else False

def better_k2(new,old):
    return True if get_scores(new)[1]<get_scores(old)[1] else False

def is_in_gini_taboo_region(coord):
    for x in taboo_gini_positions:
        if coord[0]==x[0] or coord[0]==x[1] or coord[1]==x[0] or coord[1]==x[1]:
            return True
    return False

def is_in_k2_taboo_region(coord):
    for x in taboo_k2_positions:
        if coord[0]==x[0] or coord[0]==x[1] or coord[1]==x[0] or coord[1]==x[1]:
            return True
    return False

def global_search(index_beta, prev_best_flower, prev_flower):
    prev_best = np.array(prev_best_flower)
    prev = np.array(prev_flower)
    vector = np.subtract(prev_best, prev)
    levy_vector = np.zeros(dimension)
    for i in range(dimension):
        levy_vector[i] = levy_flight(index_beta)*vector[i]
    new_flower_tmp = np.add(prev, levy_vector)
    return handle_same_SNPs(clip_agent(new_flower_tmp))

def multicriterialsort(agents,scores):
    dominated_flag = [0] * len(agents) #every agent is nondominated at start
    for i in np.arange(len(agents)):
        for j in np.arange(len(agents)):
            if (i==j) or is_the_same(agents[j]) is True:
                continue            
            if scores[j][0]<=scores[i][0] and scores[j][1]<=scores[i][1]:
                if scores[j][0]<scores[i][0] or scores[j][1]<scores[i][1]:
                    dominated_flag[i] = 1
    non_dominated_list = []
    for i in np.arange(len(agents)):
        if dominated_flag[i]==0:
            non_dominated_list.append(i)
    return dominated_flag,non_dominated_list

def clip_agent(position):
    if position[0]<lowest_coord:
        position[0] = np.random.randint(lowest_coord, highest_coord)
    if position[1]<lowest_coord:
        position[1] = np.random.randint(lowest_coord, highest_coord)
    if position[0]>highest_coord:
        position[0] = np.random.randint(lowest_coord, highest_coord)
    if position[1]>highest_coord:
        position[1] = np.random.randint(lowest_coord, highest_coord)
    return position  

def simplex_gini(xs,xg,xb):
    #get centroid of two best flowers - STEP 2
    xc = [(x + y)/2 for x, y in zip(xg, xb)]
    #make a reflection - STEP 3, xr = xc+alfa*(xc-xs)
    xr = [x+reflection_coef*(x-y) for x, y in zip(xc, xs)]
    xr = handle_same_SNPs(clip_agent(xr))
    #if xr is worse than xg, do expansion - step4 xe = xc+gama*(xr-xc)
    if better_gini(xg,xr):
        xe = [x+extension_coef*(y-x) for x, y in zip(xc, xr)]
        xe = handle_same_SNPs(clip_agent(xe))
        #if xe is worse than xg, replace xs with xr, otherwise with xe
        if better_gini(xg,xe): 
            xs = xr
        else:
            xs = xe
    
    #if xr is better than xs, do COMPRESSION xc+beta*(xs-xc) 
    if not better_gini(xs,xr):
        xt = [x+compression_coef*(y-x) for x,y in zip(xc,xs)] #COMPRESSION xc+beta*(xs-xc) 
        xt = handle_same_SNPs(clip_agent(xt))
        #if xs is better than xs
        if better_gini(xt,xs):
            xs = xt

    if better_gini(xr,xs) and better_gini(xg,xr):
        xw = [x-shrink_coef*(y-x) for x,y in zip(xc,xs)] #SHRINK xc-beta(xs-xc)
        xw = handle_same_SNPs(clip_agent(xw))
        #if xs is better than xw, than replace xs with xw, otherwise xs with xr
        if better_gini(xs,xw):
            xs = xw
        else:
            xs = xr
    return xs

def simplex_k2(xs,xg,xb):
    #get centroid of two best flowers - STEP 2
    xc = [(x + y)/2 for x, y in zip(xg, xb)]
    #make a reflection - STEP 3, xr = xc+alfa*(xc-xs)
    xr = [x+reflection_coef*(x-y) for x, y in zip(xc, xs)]
    xr = handle_same_SNPs(clip_agent(xr))
    #if xr is worse than xg, do expansion - step4 xe = xc+gama*(xr-xc)
    if better_k2(xg,xr):
        xe = [x+extension_coef*(y-x) for x, y in zip(xc, xr)]
        xe = handle_same_SNPs(clip_agent(xe))
        #if xe is worse than xg, replace xs with xr, otherwise with xe
        if better_k2(xg,xe): 
            xs = xr
        else:
            xs = xe
    
    #if xr is better than xs, do COMPRESSION xc+beta*(xs-xc) 
    if not better_k2(xs,xr):
        xt = [x+compression_coef*(y-x) for x,y in zip(xc,xs)] #COMPRESSION xc+beta*(xs-xc) 
        xt = handle_same_SNPs(clip_agent(xt))
        #if xs is better than xs
        if better_k2(xt,xs):
            xs = xt

    if better_k2(xr,xs) and better_k2(xg,xr):
        xw = [x-shrink_coef*(y-x) for x,y in zip(xc,xs)] #SHRINK xc-beta(xs-xc)
        xw = handle_same_SNPs(clip_agent(xw))
        #if xs is better than xw, than replace xs with xw, otherwise xs with xr
        if better_k2(xs,xw):
            xs = xw
        else:
            xs = xr
    return xs

class FlowerPollination:
    def __init__(self, number_of_SNPs,population_size):
        self.position_gini = generate_population(gini_population,dimension)
        self.position_k2 = generate_population(k2_population,dimension)
        self.unchanged_gini = np.zeros(gini_population)
        self.unchanged_k2 = np.zeros(k2_population)
        self.scores_gini = []
        self.scores_k2 = []
    def taboo_gini_generate(self):
        new_position = np.random.randint(lowest_coord, highest_coord,dimension)
        while is_in_gini_taboo_region(new_position) or is_the_same(new_position):
            new_position = np.random.randint(lowest_coord, highest_coord,dimension)
        return standardize_coord(new_position)
    
    def taboo_k2_generate(self):
        new_position = np.random.randint(lowest_coord, highest_coord,dimension)
        while is_in_k2_taboo_region(new_position) or is_the_same(new_position):
            new_position = np.random.randint(lowest_coord, highest_coord,dimension)
        return standardize_coord(new_position)
    
    def evaluate_scores(self):                
        self.scores_gini = [ get_scores(n)[0] for n in self.position_gini ]
        self.scores_k2 = [ get_scores(n)[1] for n in self.position_k2 ]
    
    def update(self,best_gini,best_k2,curr_iter,secondbest_gini,secondbest_k2):  
        new_positions_gini = self.position_gini.copy()
        new_positions_k2 = self.position_k2.copy()
            
        # get the value for switching between global and local search
        switch_prob = get_prob_switch(init_prob, iteration_size, curr_iter)
        probs = np.random.uniform(0, 1,(gini_population))
        do_global = probs<switch_prob
            
        for j in np.arange(gini_population):
            if do_global[j]:
                new_positions_gini[j] = global_search(index_beta, self.position_gini[best_gini],new_positions_gini[j])   
            else:
                #DO LOCAL SEARCH
                rand_indices = np.random.choice(gini_population, 2,replace=False)
                agent1, agent2 = self.position_gini[rand_indices[0]],self.position_gini[rand_indices[1]]
                new_positions_gini[j] = local_search_elite(self.position_gini[j], agent1, agent2,self.position_gini[best_gini])
                
                #DO CROSSOVER
                ec_rand = np.random.uniform(0, 1)
                if ec_rand<elite_crossover:
                    new_ls1,new_ls2 = new_positions_gini[j].copy(),new_positions_gini[j].copy()
                    rand_dim = random.randint(0,dimension-1)
                    selected = self.position_gini[np.random.choice(gini_population)]
                    new_ls1[0] = selected[rand_dim]
                    new_ls2[1] = selected[rand_dim]
                    new_ls1 = handle_same_SNPs(new_ls1)
                    new_ls2 = handle_same_SNPs(new_ls2)

                    if better_gini(new_ls1,new_ls2):
                        new_positions_gini[j] = new_ls1
                    else:
                        new_positions_gini[j] = new_ls2

            if is_in_gini_taboo_region(new_positions_gini[j]):
                new_positions_gini[j] = self.taboo_gini_generate()
            
            if better_gini(new_positions_gini[j],self.position_gini[j]):
                self.position_gini[j] = new_positions_gini[j]
                self.unchanged_gini[j] = 0
            else:
                self.unchanged_gini[j] += 1
                if self.unchanged_gini[j]>bee_limit:
                    self.unchanged_gini[j] = 0
                    self.position_gini[j] = simplex_gini(self.position_gini[j],self.position_gini[best_gini],secondbest_gini)
            
        probs = np.random.uniform(0, 1,(k2_population))
        do_global = probs<switch_prob
        for j in np.arange(k2_population):
            if do_global[j]:
                new_positions_k2[j] = global_search(index_beta, self.position_k2[best_k2],new_positions_k2[j])
            else:
                #DO LOCAL SEARCH                        
                rand_indices = np.random.choice(k2_population, 2,replace=False)
                agent1, agent2 = self.position_k2[rand_indices[0]],self.position_k2[rand_indices[1]]
                new_positions_k2[j] = local_search_elite(self.position_k2[j], agent1, agent2,self.position_k2[best_k2])
                        
                #DO CROSSOVER
                ec_rand = np.random.uniform(0, 1)
                if ec_rand<elite_crossover:
                    new_ls1,new_ls2 = new_positions_k2[j].copy(),new_positions_k2[j].copy()
                    rand_dim = random.randint(0,dimension-1)
                    selected = self.position_k2[np.random.choice(k2_population)]
                    new_ls1[0] = selected[rand_dim]
                    new_ls2[1] = selected[rand_dim]
                    new_ls1 = handle_same_SNPs(new_ls1)
                    new_ls2 = handle_same_SNPs(new_ls2)
         
                    if better_k2(new_ls1,new_ls2):
                        new_positions_k2[j] = new_ls1
                    else:
                        new_positions_k2[j] = new_ls2
                
            if is_in_k2_taboo_region(new_positions_k2[j]):
                new_positions_k2[j] = self.taboo_k2_generate()
            
            if better_k2(new_positions_k2[j],self.position_k2[j]):
                self.position_k2[j] = new_positions_k2[j]
                self.unchanged_k2[j] = 0
            else:
                self.unchanged_k2[j] += 1
                if self.unchanged_k2[j]>bee_limit:
                    self.unchanged_k2[j] = 0
                    self.position_k2[j] = simplex_k2(self.position_k2[j],self.position_k2[best_k2],secondbest_k2)
        
    def search(self):
        previous_gini, previous_k2 = [-1,-1], [-1,-1]
        prev_gini_count, prev_k2_count = 0, 0
                
        for i in np.arange(iteration_size):
            self.evaluate_scores()
            
            print("Iteration "+str(i+1),end="\r")
            #FIND THE BEST SOLUTION
            best_gini= np.argmin(self.scores_gini)
            best_gini_score = self.scores_gini[best_gini]
            best_k2 = np.argmin(self.scores_k2)
            best_k2_score = self.scores_k2[best_k2]
            
            secondbest_gini_idx = np.argpartition(self.scores_gini, 1)[1] #vybere prvych dvoch
            secondbest_gini = self.position_gini[secondbest_gini_idx]
            
            secondbest_k2_idx = np.argpartition(self.scores_k2, 1)[1] #vybere prvych dvoch
            secondbest_k2 = self.position_k2[secondbest_k2_idx]
            
            pos_gini = self.position_gini[best_gini]
            pos_gini_r = round_coord(pos_gini)
            if pos_gini_r == previous_gini:
                prev_gini_count += 1
                if prev_gini_count>required_iterations_for_taboo:
                    allSNP = []
                    for y,agent in enumerate(self.position_gini):
                        inside_position = round_coord(self.position_gini[y])
                        if is_in_taboo_coord(inside_position,pos_gini_r):
                            if inside_position[0] not in pos_gini_r:
                                allSNP.append(inside_position[0])
                            if inside_position[1] not in pos_gini_r:
                                allSNP.append(inside_position[1])
                            self.position_gini[y] = self.taboo_gini_generate()
                        else:
                            allSNP.append(inside_position[0])
                            allSNP.append(inside_position[1])
                            
                    allSNP = set(allSNP)
                    approx = best_gini_score*zeta_radius
                    for j in allSNP:
                        combined_SNP = standardize_coord([pos_gini_r[0],j])
                        if get_scores(combined_SNP)[0]<approx:
                            best_solutions.append(combined_SNP)
                        combined_SNP = standardize_coord([pos_gini_r[1],j])
                        if get_scores(combined_SNP)[0]<approx:
                            best_solutions.append(combined_SNP)
                    
                    taboo_gini_positions.append(pos_gini_r)
            else:
                prev_gini_count = 0
                previous_gini = pos_gini_r
            
            pos_k2 = self.position_k2[best_k2]
            pos_k2_r = round_coord(pos_k2)
            if pos_k2_r == previous_k2:
                prev_k2_count += 1
                if prev_k2_count>required_iterations_for_taboo:
                    allSNP = []
                    for y,agent in enumerate(self.position_k2):                
                        inside_position = round_coord(self.position_k2[y])
                        if is_in_taboo_coord(inside_position,pos_k2_r):
                            if inside_position[0] not in pos_k2_r:
                                allSNP.append(inside_position[0])
                            if inside_position[1] not in pos_k2_r:
                                allSNP.append(inside_position[1])
                            self.position_k2[y] = self.taboo_k2_generate()
                        else:
                            allSNP.append(inside_position[0])
                            allSNP.append(inside_position[1])
                      
                    allSNP = set(allSNP)
                    approx = best_k2_score*zeta_radius
                    for j in allSNP:
                        combined_SNP = standardize_coord([pos_k2_r[0],j])
                        if get_scores(combined_SNP)[1]<approx:
                            best_solutions.append(combined_SNP)
                        combined_SNP = standardize_coord([pos_k2_r[1],j])
                        if get_scores(combined_SNP)[1]<approx:
                            best_solutions.append(combined_SNP)
                            
                    taboo_k2_positions.append(pos_k2_r)
            else:
                prev_k2_count = 0  
                previous_k2 = pos_k2_r
                
            if pos_gini_r not in best_solutions:
                best_solutions.append(pos_gini_r)
            if pos_k2_r not in best_solutions:
                best_solutions.append(pos_k2_r)    
                
            self.update(best_gini,best_k2,i,secondbest_gini,secondbest_k2)  
            
with open(output_file, "w") as f:
    print("Running with the following params: ",file=f)
    print("Main path: "+path,file=f)
    print("Gini population size: "+str(gini_population),file=f)
    print("K2 population size: "+str(k2_population),file=f)
    print("number of iterations: "+str(iteration_size),file=f)
    print("Initial switch prob: "+str(init_prob),file=f)
    print("Index beta: "+str(index_beta),file=f)
    print("Iterations to ban"+str(required_iterations_for_taboo),file=f)
    print("Zeta radius: "+str(zeta_radius),file=f)
    print("*"*70,file=f)

df = pd.read_csv(path,delimiter=',')
num_of_samples = len(df)
df_without_class = df.drop('Class', axis=1)
number_of_SNPs = len(df.index)
w, h = len(df_without_class.columns), len(df_without_class.columns)
k2score_computed = [[-1 for x in range(w)] for y in range(h)] 
gini_computed = [[-1 for x in range(w)] for y in range(h)]
precomputed_p_values = [[-1 for x in range(w)] for y in range(h)]
lowest_coord = -0.49 #global parameters
highest_coord = h-0.51 #global parameters
total_combinations = (h*(h-1))/2
bonferonni = alpha_value/total_combinations

print("Detecting SNP epistasis in"+path,flush=True)

taboo_gini_positions = []
taboo_k2_positions = []
best_solutions = []
b = FlowerPollination(int(round(highest_coord)),population_size)
b.search()            
final_scores = [get_scores(j) for j in best_solutions]
flags,nondominatedsolutions = multicriterialsort(best_solutions,final_scores)    
gini_final_scores = [get_scores(j)[0] for j in best_solutions]
k2_final_scores = [get_scores(j)[1] for j in best_solutions]   
max_report = min(5,len(gini_final_scores)-1)
bestginis = np.argpartition(gini_final_scores, max_report)[:max_report+1]
max_report = min(5,len(k2_final_scores)-1)
bestk2s = np.argpartition(k2_final_scores, max_report)[:max_report+1]
top_SNPs = []
for j in bestginis:
    if best_solutions[j] not in top_SNPs:
        top_SNPs.append(best_solutions[j])
for j in bestk2s:
    if best_solutions[j] not in top_SNPs:
        top_SNPs.append(best_solutions[j]) 

for j in nondominatedsolutions:
    if best_solutions[j] not in top_SNPs:
        top_SNPs.append(best_solutions[j]) 

if not taboo_gini_positions:
    for j in taboo_gini_positions:
        if taboo_gini_positions[j] not in top_SNPs:
            top_SNPs.append(taboo_gini_positions[j])

if not taboo_k2_positions:
    for j in taboo_k2_positions:
        if taboo_k2_positions[j] not in top_SNPs:
            top_SNPs.append(taboo_k2_positions[j])

only_SNP = []
for j in top_SNPs:
    only_SNP.append(j[0])
    only_SNP.append(j[1])
only_SNP = set(only_SNP)
combs = combinations(only_SNP, 2)
comb_SNPs = [standardize_coord(round_coord(j)) for j in combs]
p_values = [ get_p_value(j) for j in comb_SNPs ]
                
max_report = min(10,len(p_values)-1)
final_best_ranks = np.argpartition(p_values, max_report)[:max_report+1]
results = [comb_SNPs[j] for j in final_best_ranks]
realSNPs = []
for j in results:
    realSNP = []
    realSNP.append(df.columns[j[0]])
    realSNP.append(df.columns[j[1]])
    realSNPs.append(realSNP)

print("Results of the 1st stage of FPepi algorithm",realSNPs)    
p_values = [ get_p_value(j) for j in results ]
        
print("Final results of the FPepi algorithm",realSNPs)    
final_prisne = [results[index] for index,j in enumerate(p_values) if j<bonferonni ]         
if not final_prisne:
    print("No SNP combination passed the G-test",flush=True)
else:
    print("The candidate set after the G-test",flush=True)
    for j in final_prisne:
        print(df.columns[j[0]],df.columns[j[1]],"p_value:",get_p_value(j),"Gini score:",get_scores(j)[0],"K2 score:",get_scores(j)[1])
                                       
#print  to file
with open(output_file, "a") as f:
    print("Results of FPepi for "+ path,file=f)
    print("Results of the 1st stage of Fpepi algorithm:",file=f)
    print(realSNPs,file=f)
            
    print("Final results after the G-test with p_value",bonferonni,file=f)
    if not final_prisne:
        print("No SNP combination passed the G-test",file=f)
    else:
        for j in final_prisne:
            realSNP = []
            realSNP.append(df.columns[j[0]])
            realSNP.append(df.columns[j[1]])
            print("    ",realSNP,"p_value:",get_p_value(j),file=f)