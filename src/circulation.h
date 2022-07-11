#ifndef CIRCULATION_H
#define CIRCULATION_H

#include "defines.h"
#include "graph.h"
#include "tree.h"

struct nvalue_t;
struct evalue_t;

typedef node_t<nvalue_t, evalue_t> vertex_t;
typedef edge_t<nvalue_t, evalue_t> arc_t;

typedef node_t<vertex_t *, arc_t *> vertexp_t;
typedef edge_t<vertex_t *, arc_t *> arcp_t;

TAILQ_HEAD(vertexhead, vertex_t);


struct dist_t {
	dist_t() : dist(0), depth(0), bias(0) {}
	weight_t dist;  // dual distance
	uint32_t depth; // prefer shallow trees (also, needed for rr algorithm to work)
	double bias;    // in case of a tie, go for the most biased node

	bool
	operator < (const dist_t & a)
	const 
	{
		return (dist < a.dist) || (dist == a.dist && depth < a.depth) || (dist == a.dist && depth == a.depth && bias < a.bias);
	}

	bool operator > (const dist_t & a) const { return a < *this;}

	bool operator == (const dist_t & a) const { return dist == a.dist && depth == a.depth && bias == a.bias;}

};

struct nvalue_t {
	nvalue_t()
	{
		circ.bias = 0; circ.dual = 0; circ.flow = 0;
		rr.isaffected = rr.inworkset = rr.instack = false;
	}

	// if this is in constructor, it will be invalidated during copying.
	void
	init(vertex_t *v)
	{
		contract.master = v;
		TAILQ_INIT(&contract.slaves);
	}

	uint32_t label;

	struct {
		vertexhead slaves;
		TAILQ_ENTRY(vertex_t) slave;
		vertex_t *master;

		bool instack;
		TAILQ_ENTRY(vertex_t) stack;
	} contract;

	struct {
		double bias;
		double flow;
		weight_t dual;
	} circ;

	struct {
		TAILQ_ENTRY(vertex_t) imbalance, affected, workset;
		
		dist_t distance;

		bool isaffected;
		bool completed;
		bool inworkset;
		bool instack;

		vertexp_t *node;

		RB_ENTRY(vertex_t) stack;


	} rr; // RR algorithm for dynamic shortest tree
};

RB_HEAD(vertextree, vertex_t);


struct evalue_t {
	evalue_t() {circ.flow = circ.cost = 0;}
	struct {
		double flow;
		weight_t cost;
	} circ;

	struct {
		arcp_t *f, *b;
	} rr;
};


typedef graph<nvalue_t, evalue_t> network;
typedef graph<vertex_t *, arc_t *> residual;

typedef std::vector<network *> networkvector; 

struct compnode_t;
struct compedge_t;

typedef node_t<compnode_t, compedge_t> cvertex_t;
typedef edge_t<compnode_t, compedge_t> carc_t;

TAILQ_HEAD(cvertexhead, cvertex_t);

struct compnode_t {
	compnode_t() : lowlink(0), index(0), instack(false), compid(0), compind(0) {}

	uint32_t label;

	uint32_t lowlink;
	uint32_t index;
	bool instack;

	uint32_t dual;

	uint32_t compid;
	uint32_t compind;

	TAILQ_ENTRY(cvertex_t) stack, root, comp;
	cvertexhead comps;
};

struct compedge_t {
	double cost;
	TAILQ_ENTRY(carc_t) stack;
};

typedef graph<compnode_t, compedge_t> compgraph;

compgraph * read(FILE *f, bool weighted, bool self);
network *read(FILE *f, bool weighted, bool self, uint32_t & ncnt, uint32_t & ecnt, uint32_t limit);

void balance(network & g, residual & r, double delta);
void canonize(network & g, residual & r, vertex_t *canon);
void canonize(network & g, residual & r);
void init_residual(network & g, residual & r);

vertex_t * computeblocks(network & n, double threshold);
void batchcontract(network & g, vertex_t *first, network::edgevector & edges);
void unroll_dual(network & g);
void unroll_master(network & g);


uint32_t ccomp(compgraph & g, cvertexhead *roots);
network * create_network(compgraph & g, cvertex_t *root, uint32_t & ncnt, uint32_t & ecnt, uint32_t limit);
network * join(networkvector & comps, compgraph & cg);
void migrate_dual(compgraph & cg, networkvector & comps, network & j);

#endif
