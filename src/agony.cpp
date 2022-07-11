#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "defines.h"
#include "graph.h"
#include "circulation.h"



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

		printf("delta: %f %f %f\n", delta, d1, d2);


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
			printf("FINAL: %d %d %.2f %.2f %d\n", v->id, v->v.circ.dual, v->v.circ.bias,  v->v.circ.flow, v->v.contract.master);
		}*/
	}

	unroll_dual(g);

	/*
	for (uint32_t i = 0; i < g.nodebudget(); i++) {
		vertex_t *v = g.get(i);
		printf("DONE: %d %d %f\n", v->id, v->v.circ.dual, v->v.circ.bias - v->v.circ.flow);
	}*/
}

void
outputrank(FILE *f, network & g, uint32_t ncnt)
{
	for (uint32_t i = 0; i < ncnt; i++) {
		vertex_t *v = g.get(i);
		fprintf(f, "%d %d\n", v->v.label, v->v.circ.dual);
	}
}

void
outputrank(FILE *f, compgraph & g, uint32_t ncnt)
{
	for (uint32_t i = 0; i < ncnt; i++) {
		cvertex_t *v = g.get(i);
		fprintf(f, "%d %d\n", v->v.label, v->v.dual);
	}
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



int
main(int argc, char **argv)
{
	static struct option longopts[] = {
		{"out",             required_argument,  NULL, 'o'},
		{"in",              required_argument,  NULL, 'i'},
		{"groups",          required_argument,   NULL, 'k'},
		{"weighted",        no_argument,        NULL, 'w'},
		{"self",            no_argument,        NULL, 's'},
		{"noscc",           no_argument,        NULL, 'C'},
		{"help",            no_argument,        NULL, 'h'},
		{ NULL,             0,                  NULL,  0 }
	};

	char *inname = NULL;
	char *outname = NULL;
	bool weighted = false;
	bool self = false;
	bool scc = true;
	uint32_t groupcnt = 0;


	int ch;
	while ((ch = getopt_long(argc, argv, "whk:o:i:sC", longopts, NULL)) != -1) {
		switch (ch) {
			case 'h':
				printf("Usage: %s -i <input file> -o <output file> [-k <number>] [-whsC]\n", argv[0]);
				printf("  -h    print this help\n");
				printf("  -i    input file\n");
				printf("  -o    output file\n");
				printf("  -k    number of groups\n");
				printf("  -w    input file has weights\n");
				printf("  -s    allow self loops\n");
				printf("  -C    disable strongly connected components\n");
				return 0;
				break;
			case 'i':
				inname = optarg;
				break;
			case 'o':
				outname = optarg;
				break;
			case 'w':
				weighted = true;
				break;
			case 'k':
				groupcnt = atoi(optarg);
				break;
			case 's':
				self = true;
				break;
			case 'C':
				scc = false;
				break;
		}

	}
	if (inname == NULL) {
		printf("Missing input file\n");
		return 1;
	}

	uint32_t ncnt; 
	uint32_t ecnt; 
	


	if (!scc || groupcnt > 0) {
		FILE *f = fopen(inname, "r");
		network *g = read(f, weighted, self, ncnt, ecnt, groupcnt);
		fclose(f);

		residual r(g->nodecnt(), 2*g->edgebudget());
		init_residual(*g, r);
		vertex_t *canon = g->get(g->nodecnt() - 2);

		computeagony(*g, r, canon);
		printf("Agony: %f \n", score(*g, ncnt, ecnt));

		f = stdout;
		if (outname) f = fopen(outname, "w");
		outputrank(f, *g, ncnt);
		if (outname) fclose(f);

		delete g;
	}
	else {

		FILE *f = fopen(inname, "r");
		compgraph *cg = read(f, weighted, self);
		fclose(f);

		cvertexhead roots;
		uint32_t compcnt = ccomp(*cg, &roots);
		
		printf("components: %d\n", compcnt);

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
		printf("Agony: %f \n", sc);

		network *joined = join(comps, *cg);
		residual r(joined->nodecnt(), 2*joined->edgebudget());
		init_residual(*joined, r);

		canonize(*joined, r);
		migrate_dual(*cg, comps, *joined);

		f = stdout;
		if (outname) f = fopen(outname, "w");
		outputrank(f, *cg, cg->nodebudget());
		if (outname) fclose(f);

		for (uint32_t i = 0; i < comps.size(); i++)
			delete comps[i];

		delete joined;
		delete cg;
	}
	



	return 0;
}
