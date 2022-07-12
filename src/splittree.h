#ifndef SPLITTREE_H
#define SPLITTREE_H

#include "network.h"
#include <list>


struct splitvertex;


struct splitvertex {

	splitvertex() :
		back(0), inback(0), outback(0), diffpos(0), diffneg(0),
			positive(0), negative(0), positive_zero(0), negative_zero(0),
		left(0), right(0), gain(0),
		budget(0), rank(0) {}

	weight_t back;
	weight_t inback;
	weight_t outback;
	weight_t diffpos;
	weight_t diffneg;

	vertexhead *positive;
	vertexhead *negative;
	vertexhead *positive_zero;
	vertexhead *negative_zero;

	struct splitvertex *left;
	struct splitvertex *right;

	weight_t gain;

	// dynamic program
	weightvector opt;
	uintvector ind;

	uint32_t budget;
	uint32_t rank;

};

typedef std::list<splitvertex *> svlist;

splitvertex * init(vertexhead *h, svlist & tree);


void split(splitvertex *a, network & g, svlist & tree);

void splitleft(splitvertex *a, network & g, svlist & tree);
void splitright(splitvertex *a, network & g, svlist & tree);

splitvertex * init(network & g, svlist & tree);
void free(svlist & tree);

void init_contract(svlist & tree, uint32_t k);
void compute_contract(svlist & tree);
void set_restricted_budget(svlist & tree, uint32_t k);

uint32_t set_unrestricted_budget(svlist & tree);

void set_ranks(svlist & tree, uint32_t base);

#endif
