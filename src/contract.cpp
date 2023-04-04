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

#include "defines.h"
#include "graph.h"
#include "circulation.h"

#include <stdio.h>
#include <algorithm>


vertex_t *
computeblocks(network & n, double threshold)
{
	vertexhead notvisited;
	TAILQ_INIT(&notvisited);

	for (vertex_t *v = n.get_first(); v; v = v->next()) {
		v->v.contract.instack = false;
		TAILQ_INSERT_TAIL(&notvisited, v, v.contract.stack);
	}

	vertex_t *first = 0;

	// Do bfs to find connected components with edge flow >= threshold
	while (!TAILQ_EMPTY(&notvisited)) {
		vertex_t *r = TAILQ_FIRST(&notvisited);
		TAILQ_REMOVE(&notvisited, r, v.contract.stack);

		// check if the node is isolated
		bool isolated = true;
		for (arc_t *e = r->first_out(); e && isolated; e = e->next_out())
			if (e->v.circ.flow >= threshold) isolated = false;
		for (arc_t *e = r->first_in(); e && isolated; e = e->next_in())
			if (e->v.circ.flow >= threshold) isolated = false;
		
		if (isolated) continue;


		vertexhead stack; 
		TAILQ_INIT(&stack);
		TAILQ_INSERT_TAIL(&stack, r, v.contract.stack);
		r->v.contract.instack = true;

		vertex_t *master = n.addnode();
		master->v.init(master);
		if (!first) first = master;

		//printf("root: %d, master %d\n", r->id, master->id);

		while (!TAILQ_EMPTY(&stack)) {
			vertex_t *v = TAILQ_FIRST(&stack);
			//printf("vertex: %d\n", v->id);
			TAILQ_REMOVE(&stack, v, v.contract.stack);

			TAILQ_INSERT_TAIL(&master->v.contract.slaves, v, v.contract.slave);
			v->v.contract.master = master;

			for (arc_t *e = v->first_out(); e; e = e->next_out()) {
				vertex_t *w = e->to;
				if (e->v.circ.flow >= threshold && !w->v.contract.instack) {
					w->v.contract.instack = true;
					TAILQ_REMOVE(&notvisited, w, v.contract.stack);
					TAILQ_INSERT_TAIL(&stack, w, v.contract.stack);
				}
			}

			for (arc_t *e = v->first_in(); e; e = e->next_in()) {
				vertex_t *w = e->from;
				if (e->v.circ.flow >= threshold && !w->v.contract.instack) {
					w->v.contract.instack = true;
					TAILQ_REMOVE(&notvisited, w, v.contract.stack);
					TAILQ_INSERT_TAIL(&stack, w, v.contract.stack);
				}
			}
		}
	}
	return first;
}

bool
edgecmp(arc_t *a, arc_t *b)
{
	uint32_t a1 = a->from->v.contract.master->id;
	uint32_t b1 = b->from->v.contract.master->id;
	uint32_t a2 = a->to->v.contract.master->id;
	uint32_t b2 = b->to->v.contract.master->id;

	if (a1 != b1) return a1 < b1;
	if (a2 != b2) return a2 < b2;
	return a->v.circ.cost < b->v.circ.cost;
}

void
batchcontract(network & g, vertex_t *first, network::edgevector & edges)
{
	uint32_t ecnt = 0;

	// collect all edges, set dual and bias for new nodes
	for (vertex_t *m = first; m; m = m->next()){
		vertex_t *s;

		m->v.circ.dual = TAILQ_FIRST(&m->v.contract.slaves)->v.circ.dual; 

		TAILQ_FOREACH(s, &m->v.contract.slaves, v.contract.slave) {

			m->v.circ.bias += s->v.circ.bias; // compute new bias as the sum of slaves' bias
			m->v.circ.flow += s->v.circ.flow; // compute new flow as the sum of slaves' flow

			for (arc_t *e = s->first_out(); e; e = e->next_out()) {
				edges[ecnt++] = e;
			}

			for (arc_t *e = s->first_in(); e; e = e->next_in()) {
				vertex_t *w = e->from;
				if (w->v.contract.master->id < first->id) // edges incoming from singletons to new blocks
					edges[ecnt++] = e;
			}
		}
	}

	// sort edges
	std::sort(edges.begin(), edges.begin() + ecnt, edgecmp);

	// reassign edges
	vertex_t *prevfrom = 0, *prevto = 0;
	arc_t *prevedge;
	for (uint32_t i = 0; i < ecnt; i++) {
		arc_t *e = edges[i];
		vertex_t *from = e->from->v.contract.master;
		vertex_t *to = e->to->v.contract.master;
		if (from == to) {
			g.unbindedge(e); // delete loop edges
		}
		else { // if (from != prevfrom || to != prevto) {
			// new cost (note that it may be negative but the duals will make it positive)
			e->v.circ.cost = e->v.circ.cost - (e->from->v.circ.dual - from->v.circ.dual) + (e->to->v.circ.dual - to->v.circ.dual);
			// rebind edge
			g.bindedge(e, from, to);
			prevedge = e;
		}
		/*
		XXX clean this code, this used to be a bug
		we cannot just add the flow because we may have edges with different costs
		else {
			// add flow
			prevedge->v.circ.flow += e->v.circ.flow;
			g.unbindedge(e); // delete redundant edge
		}
		*/
		prevfrom = from;
		prevto = to;
	}

	// delete slave nodes
	for (vertex_t *m = first; m; m = m->next()){
		vertex_t *s;

		TAILQ_FOREACH(s, &m->v.contract.slaves, v.contract.slave) {
			s->v.circ.dual -= m->v.circ.dual; // update slaves' dual
			g.deactivatenode(s);
		}
	}
}

void
unroll_dual(network & g)
{
	for (int32_t i = g.nodebudget() - 1; i >= 0; i--) {
		vertex_t *m = g.get(i);
		vertex_t *s;
		//printf("master: %d\n", m->id);
		TAILQ_FOREACH(s, &m->v.contract.slaves, v.contract.slave) {
			//printf("slave: %d\n", s->id);
			s->v.circ.dual += m->v.circ.dual; // update slaves' dual
		}
	}
}


// flattens the master tree structure, so that the leaves (of componets) point
// out to the root. after this masters and slaves list do not match 
void
unroll_master(network & g)
{
	for (int32_t i = g.nodebudget() - 1; i >= 0; i--) {
		vertex_t *m = g.get(i);
		vertex_t *s;
		TAILQ_FOREACH(s, &m->v.contract.slaves, v.contract.slave) {
			s->v.contract.master = m->v.contract.master;
		}
	}
}
