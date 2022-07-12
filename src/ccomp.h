#ifndef CCOMP_H
#define CCOMP_H

#include "network.h"
#include "splittree.h"


struct sccnvalue_t;
struct sccevalue_t;

typedef node_t<sccnvalue_t, sccevalue_t> sccvertex_t;
typedef edge_t<sccnvalue_t, sccevalue_t> sccarc_t;


TAILQ_HEAD(sccvertexhead, sccvertex_t);
TAILQ_HEAD(sccarchead, sccarc_t);


struct sccnvalue_t {
	sccnvalue_t() {}

	svlist tree;
	vertexhead elems;

    // dynamic program
	weightvector opt;
	uintvector ind;
	uintvector cluster;

	uint32_t budget;
};


struct sccevalue_t {
	sccevalue_t() {w = 0;}
	weight_t w;
};


typedef graph<sccnvalue_t, sccevalue_t> sccnetwork;


uint32_t findscc(network & g, vertexhead *order);
sccnetwork * scctolayers(network & g, vertexhead *order, uint32_t compcnt);

void init_contract(sccnetwork & g, uint32_t k);
void compute_contract(sccnetwork & g, uint32_t k);
void set_restricted_budget(sccnetwork & g, uint32_t k);

sccnetwork * dummyscc(network & g);



#endif
