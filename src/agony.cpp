/*
 * 
 * 
 * The software for agony computation (exact algorithm) adopted
 * here was developed by Professor Nikolaj Tatti and colleagues.
 * For a detailed description, please refer to: Tatti, Nikolaj.
 * "Tiers for peers: a practical algorithm for discovering hierarchy
 * in weighted networks." Data mining and knowledge discovery
 * 31.3 (2017): 702-738.
 *
 * Copyright 2002 Niels Provos <provos@citi.umich.edu>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "defines.h"
#include "graph.h"
#include "circulation.h"
#include <Rcpp.h>
using namespace Rcpp;



void
computeagony(network & g, residual & r, vertex_t *canon)
{
	double delta = 0;

	uint32_t ncnt = g.nodecnt();

	network::edgevector edges(g.edgebudget());

	while (g.nodecnt() > 1) {
		double d1 = 0, d2 = 0;
		for (vertex_t *v = g.get_first(); v; v = v->next()) {
			d1 = std::max(d1, v->v.circ.bias - v->v.circ.flow);
			d2 = std::min(d2, v->v.circ.bias - v->v.circ.flow);
			//printf("%d %f\n", v->v.circ.dual, v->v.circ.bias);
		}
		double d = std::min(d1, -d2);
		if (d == 0) break;
		if (delta == 0) {
			int exp;
			if (frexp(d, &exp) == 0.5) exp--;
			delta = ldexp(1.0, exp);
		}
		else {
			while (delta*3/4 > d) delta /= 2;
		}

		//printf("delta: %f %f %f\n", delta, d1, d2);


		if (g.nodecnt() == 1) break;

		balance(g, r, delta);
		//printf("canon: %d\n", canon->id);
		canonize(g, r, canon);

		vertex_t *master = computeblocks(g, 3*ncnt*delta);
		if (master) {
			batchcontract(g, master, edges);
			init_residual(g, r); // update residual graph
			canon = canon->v.contract.master; // update canon, if changed
		}

		/*
		for (uint32_t i = 0; i < g.nodebudget(); i++) {
			vertex_t *v = g.get(i);
			printf("FINAL: %d %.2f %.2f\n", v->v.circ.dual, v->v.circ.bias,  v->v.circ.flow);
		}
		*/
	}

	unroll_dual(g);

	/*
	for (uint32_t i = 0; i < g.nodebudget(); i++) {
		vertex_t *v = g.get(i);
		printf("DONE: %d %f\n", v->v.circ.dual, v->v.circ.bias - v->v.circ.flow);
	}
	*/

}

void
outputrank(FILE *f, network & g, uint32_t ncnt)
{
	for (uint32_t i = 0; i < ncnt; i++) {
		vertex_t *v = g.get(i);
		//fprintf(f, "%d %d\n", v->v.label, v->v.circ.dual);
	}
}


// ytfyt
IntegerMatrix
outputrank(compgraph & g, uint32_t ncnt)
{
	IntegerMatrix output(ncnt, 2);
	for (uint32_t i = 0; i < ncnt; i++) {
		cvertex_t *v = g.get(i);
		output(i, 0) = v->v.label;
		output(i, 1) = v->v.dual;
	}
	return(output);
}

double
score(network & g, uint32_t ncnt, uint32_t ecnt)
{
	double agony = 0;
	for (uint32_t i = 0; i < ncnt + ecnt; i++) {
		vertex_t *v = g.get(i);
		agony += v->v.circ.dual*v->v.circ.bias;
	}
	double shift = 0;
	for (uint32_t i = 0; i < ncnt; i++) {
		vertex_t *v = g.get(i);
		shift += v->v.circ.bias;
	}

	return shift - agony;
}


//' Perform agony computation. The software for agony computation 
//' (exact algorithm) adopted here was developed by Professor Nikolaj
//' Tatti and colleagues.
//' Agony is a measure of hierarchy within a directed graph.
//' Given a directed graph and a ranking metric (e.g., in our
//' case the time ordering of accumulation of driver alterations
//' during tumor evolution), any arc from nodes that are higher
//' in the hierarchy (e.g., alterations that occur in later stages
//' of the tumor) to nodes that are lower in the hierarchy
//' (e.g., alterations that occur at the initiation of the tumor)
//' are not expected and they are said to be causing agony.
//'
//' For a detailed description, please refer to: Tatti, Nikolaj.
//' "Tiers for peers: a practical algorithm for discovering hierarchy
//' in weighted networks." Data mining and knowledge discovery
//' 31.3 (2017): 702-738.
//'
//' @title agony
//''
//' @examples
//' data(agonyArcs)
//' agony(agonyArcs, 12345)
//'
//' @param inmatrix Input agony matrix.
//' @param seed Input seed
//' @return Output agony matrix.
//'
//' 
//' The software for agony computation (exact algorithm) adopted
//' here was developed by Professor Nikolaj Tatti and colleagues.
//' For a detailed description, please refer to: Tatti, Nikolaj.
//' "Tiers for peers: a practical algorithm for discovering hierarchy
//' in weighted networks." Data mining and knowledge discovery
//' 31.3 (2017): 702-738.
//'
// [[Rcpp::export]]
Rcpp::IntegerMatrix agony(Rcpp::IntegerMatrix inmatrix, uint32_t seed)
{
  
  // seed can be used as seed.
  // let's decide if we want to use it or not
	//printf("seed: %d \n", seed);

  //bool weighted = false;
  //bool self = false;
  uint32_t groupcnt = 0;
  
  uint32_t ncnt; 
  uint32_t ecnt; 
  
  //const char* test_out = outname.get_cstring();
  
  IntegerMatrix m = as<IntegerMatrix>(inmatrix);
  
  compgraph *cg = read(m);

  cvertexhead roots;
  uint32_t compcnt = ccomp(*cg, &roots);
  
  double sc = 0;
  
  networkvector comps(compcnt);
  
  uint32_t ind = 0;
  cvertex_t *root;
  TAILQ_FOREACH(root, &roots, v.root) {
    comps[ind] = create_network(*cg, root, ncnt, ecnt, groupcnt);
    network *g = comps[ind];
    
    residual r(g->nodecnt(), 2*g->edgebudget());
    init_residual(*g, r);
    vertex_t *canon = g->get(g->nodecnt() - 2);
    
    computeagony(*g, r, canon);
    sc += score(*g, ncnt, ecnt);
    unroll_master(*g);
    ind++;
  }
  
  network *joined = join(comps, *cg);
  residual r(joined->nodecnt(), 2*joined->edgebudget());
  init_residual(*joined, r);
  
  canonize(*joined, r);
  migrate_dual(*cg, comps, *joined);
  
  IntegerMatrix output = outputrank(*cg, cg->nodebudget());
  
  for (uint32_t i = 0; i < comps.size(); i++)
    delete comps[i];
  
  delete joined;
  delete cg;

  return(output);
  
}
