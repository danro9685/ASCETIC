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
	//printf("outputrank:\n");
	for (uint32_t i = 0; i < ncnt; i++) {
		cvertex_t *v = g.get(i);
		output(i, 0) = v->v.label;
		output(i, 1) = v->v.dual;
		//printf("%d %d\n", v->v.label, v->v.dual);
		//fprintf(f, "%d %d\n", v->v.label, v->v.dual);
	}
	//printf("\n");
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


//' Perform agony computation. The software for agony computation (exact algorithm) adopted here was developed by Professor Nikolaj Tatti and colleagues.
//' Agony is a measure of hierarchy within a directed graph. Given a directed graph and a ranking metric (e.g., in our case the time ordering of accumulation 
//' of driver alterations during tumor evolution), any arc from nodes that are higher in the hierarchy (e.g., alterations that occur in later stages of the tumor) 
//' to nodes that are lower in the hierarchy (e.g., alterations that occur at the initiation of the tumor) are not expected and they are said to be causing agony.
//' For a detailed description, please refer to: Tatti, Nikolaj. "Tiers for peers: a practical algorithm for discovering hierarchy in weighted networks." Data mining and knowledge discovery 31.3 (2017): 702-738.
//'
//' @title agony
//' @param inmatrix Input agony matrix.
//' @param outname Output agony file.
//' @return Output agony matrix.
//'
// [[Rcpp::export]]
Rcpp::IntegerMatrix agony(Rcpp::IntegerMatrix inmatrix)
{
  
  //bool weighted = false;
  //bool self = false;
  uint32_t groupcnt = 0;
  
  uint32_t ncnt; 
  uint32_t ecnt; 
  
  //const char* test_out = outname.get_cstring();
  
  IntegerMatrix m = as<IntegerMatrix>(inmatrix);
  
  //compgraph *cg = read(m, weighted, self);
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
  
  //f = stdout;
  //FILE *f = fopen(test_out, "w");
  IntegerMatrix output = outputrank(*cg, cg->nodebudget());
  //fclose(f);
  
  for (uint32_t i = 0; i < comps.size(); i++)
    delete comps[i];
  
  delete joined;
  delete cg;

  return(output);
  
}
