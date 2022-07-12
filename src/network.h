#ifndef CIRCULATION_H
#define CIRCULATION_H

#include "defines.h"
#include "graph.h"

struct nvalue_t;
struct evalue_t;

typedef node_t<nvalue_t, evalue_t> vertex_t;
typedef edge_t<nvalue_t, evalue_t> arc_t;


TAILQ_HEAD(vertexhead, vertex_t);
TAILQ_HEAD(archead, arc_t);


struct nvalue_t {
	nvalue_t()
	{
		label = rank = deg = 0;
		flux = inback = outback = 0;
		scc.layer = 0;
		scc.index = scc.lowlink = 0;
		scc.instack = false;
	}

	uint32_t label;
	uint32_t rank;
	uint32_t deg;
	weight_t flux;
	weight_t inback;
	weight_t outback;

	weight_t diff() const {return flux + inback - outback;}

	struct {
		uint32_t compid;
		uint32_t layer;
		uint32_t index, lowlink;
		bool instack;
	} scc;

	TAILQ_ENTRY(vertex_t) group;

};


struct evalue_t {
	evalue_t() {w = 0; cross = false;}
	weight_t w;
	bool cross;

	TAILQ_ENTRY(arc_t) group;
};


typedef graph<nvalue_t, evalue_t> network;

uint32_t findscc(network & g);

network *read(FILE *f, bool weighted, bool self, uint32_t & ncnt, uint32_t & ecnt);

#endif
