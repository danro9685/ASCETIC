#include "circulation.h"
#include <limits>
#include <assert.h>
#include <stdio.h>

inline int
treecomp(const vertex_t *a, const vertex_t *b)
{
	if (a->id == b->id) return 0;
	if (a->v.rr.distance < b->v.rr.distance) return -1;
	if (a->v.rr.distance > b->v.rr.distance) return 1;
	return a->id < b->id ? -1 : 1;
}

RB_PROTOTYPE(vertextree, vertex_t, v.rr.stack, treecomp);
RB_GENERATE(vertextree, vertex_t, v.rr.stack, treecomp);

void
init_residual(network & g, residual & r)
{
	r.wipe();
	for (vertex_t *v = g.get_first(); v; v = v->next()) {
		vertexp_t *x = r.addnode();
		v->v.rr.node = x;
		x->v = v;

		for (arc_t *e = v->first_out(); e; e = e->next_out()) {
			arcp_t *x = r.addedge();
			arcp_t *y = r.addedge();
			y->v = x->v = e;
			e->v.rr.f = x;
			e->v.rr.b = y;
		}
	}
}

bool
augmentpath(residual & r, vertex_t *sink, double delta, double thresh, vertexhead *workset) 
{
	TAILQ_INIT(workset);
	
	
	vertexp_t *v = sink->v.rr.node;
	while (true) {
		//printf("augment step: %d \n", v->v->id);
		arcp_t *e = v->first_in();
		if (e == 0) break;
		vertexp_t *w = e->from;

		if (e == e->v->v.rr.f) { // forward edge
			e->v->v.circ.flow += delta;
		}
		else {
			e->v->v.circ.flow -= delta;
			if (e->v->v.circ.flow <= delta / 2) { // this test is equivalent to flow = 0 but done this way due to numerical stability
				e->v->v.circ.flow = 0; // just to be sure 
				r.unbindedge(e);
				assert(v->v->v.rr.inworkset == false);
				//printf("adding to workset: %d %d\n", v->v->id, w->v->id);
				TAILQ_INSERT_TAIL(workset, v->v, v.rr.workset);
				v->v->v.rr.inworkset = true;
			}
		}
		v = w;
	}

	vertex_t *root = v->v;

	assert(sink != root);

	//printf("augment: %d %d (%f %f)\n", root->id, sink->id, root->v.circ.bias - root->v.circ.flow, sink->v.circ.bias - sink->v.circ.flow);

	//if (root->v.circ.bias - root->v.circ.flow < thresh) printf("sink: %d", sink->id);
	
	assert(root->v.circ.bias - root->v.circ.flow >= thresh); 

	sink->v.circ.flow -= delta;
	root->v.circ.flow += delta;


	// no longer a root, should delete the whole tree.
	if (root->v.circ.bias - root->v.circ.flow < thresh) {
		assert(root->v.rr.inworkset == false);
		TAILQ_INSERT_TAIL(workset, root, v.rr.workset);
		root->v.rr.inworkset = true;
		return true;
	}

	return false;
}

void
processworkset(residual & r, vertexhead *workset, vertexhead *affected)
{
	TAILQ_INIT(affected);

	while (!TAILQ_EMPTY(workset)) {
		vertex_t *v = TAILQ_FIRST(workset);
		TAILQ_REMOVE(workset, v, v.rr.workset);
		assert(v->v.rr.inworkset);
		v->v.rr.inworkset = false;

		//printf("workset %d\n", v->id);

		if (v->v.rr.node->first_in() == 0) { // in_degree = 0 
			v->v.rr.isaffected = true;
			TAILQ_INSERT_TAIL(affected, v, v.rr.affected);

			arcp_t *e;
			while ((e = v->v.rr.node->first_out()) != 0) {
				vertex_t *w = e->to->v;
				//printf("adding %d\n", w->id);
				if (!w->v.rr.inworkset) {
					TAILQ_INSERT_TAIL(workset, w, v.rr.workset);
					w->v.rr.inworkset = true;
				}
				r.unbindedge(e);
			}
			assert(v->v.rr.node->first_out() == 0);
		}
	}

	vertex_t *w;
	TAILQ_FOREACH(w, affected, v.rr.affected) {
		assert(w->v.rr.node->first_out() == 0);
	}
	//printf("done\n");
}

void
updatedist(residual & r, vertextree *stack, vertex_t *v, vertex_t *w, arc_t *e)
{
	if (w->v.rr.completed) return;
	assert(!v->v.rr.isaffected);
	assert(v->v.rr.node->v == v);

	dist_t cand = v->v.rr.distance;
	assert(e->v.circ.cost + e->to->v.circ.dual - e->from->v.circ.dual >= 0);
	//printf("%d %d: %d %d %d (%d + %d)\n", v->id, w->id, v->v.circ.dual, w->v.circ.dual, e->v.circ.cost, cand.dist, e->v.circ.cost + e->to->v.circ.dual - e->from->v.circ.dual);
	cand.dist += e->v.circ.cost + e->to->v.circ.dual - e->from->v.circ.dual;
	cand.depth++;
	cand.bias = w->v.circ.bias - w->v.circ.flow;

	arcp_t *re = e->from == v ? e->v.rr.f : e->v.rr.b;

	if (cand < w->v.rr.distance) {
		if (w->v.rr.instack) RB_REMOVE(vertextree, stack, w);
		w->v.rr.distance = cand;
		RB_INSERT(vertextree, stack, w);
		w->v.rr.instack = true;
		r.unbind_incoming(w->v.rr.node);
		r.bindedge(re, v->v.rr.node, w->v.rr.node);
	}
	else if (cand == w->v.rr.distance) {
		r.bindedge(re, v->v.rr.node, w->v.rr.node);
	}
}

void
processaffected(residual & r, vertexhead *affected, vertextree *stack)
{
	vertex_t *v;
	TAILQ_FOREACH(v, affected, v.rr.affected) {
		assert(v->v.rr.node->first_out() == 0);
	}
	TAILQ_FOREACH(v, affected, v.rr.affected) {
		//printf("affected %d\n", v->id);
		assert(v->v.rr.isaffected);
		v->v.rr.distance.dist = std::numeric_limits<weight_t>::max();
		v->v.rr.completed = false;
		assert(v->v.rr.node->first_out() == 0);
		for (arc_t *e = v->first_out(); e; e = e->next_out()) {
			if (e->v.circ.flow == 0) continue;
			vertex_t *w = e->to;
			//printf("parent %d\n", w->id);
			if (w->v.rr.isaffected) continue;
			updatedist(r, stack, w, v, e);
		}

		for (arc_t *e = v->first_in(); e; e = e->next_in()) {
			vertex_t *w = e->from;
			//printf("parent %d\n", w->id);
			if (w->v.rr.isaffected) continue;
			updatedist(r, stack, w, v, e);
		}
	}

	TAILQ_FOREACH(v, affected, v.rr.affected) {
		v->v.rr.isaffected = false;
	}

}

void
processstack(residual & r, vertextree *stack)
{
	while (!RB_EMPTY(stack)) {
		vertex_t *v = RB_MIN(vertextree, stack);
		RB_REMOVE(vertextree, stack, v);
		v->v.rr.instack = false;
		assert(v->v.rr.completed == false);
		v->v.rr.completed = true;
		assert(v->v.rr.node->first_out() == 0);
		//printf("stack: %d %d\n", v->id, v->v.rr.distance.dist);

		for (arc_t *e = v->first_out(); e; e = e->next_out()) {
			vertex_t *w = e->to;
			updatedist(r, stack, v, w, e);
		}

		for (arc_t *e = v->first_in(); e; e = e->next_in()) {
			if (e->v.circ.flow == 0) continue;
			vertex_t *w = e->from;
			updatedist(r, stack, v, w, e);
		}

		// update duals
		v->v.circ.dual -= v->v.rr.distance.dist;
		v->v.rr.distance.dist = 0;
	}

	//printf("done\n");

}


void
balance(network & g, residual & r, double delta)
{
	double thresh = delta * 3 / 4;

	vertexhead sources, sinks;
	TAILQ_INIT(&sources);
	TAILQ_INIT(&sinks);
	uint32_t source_cnt = 0;

	vertextree stack;
	RB_INIT(&stack);

	for (vertex_t *v = g.get_first(); v; v = v->next()) {
		v->v.rr.instack = false;
		v->v.rr.inworkset = false;
		v->v.rr.isaffected = false;
		v->v.rr.completed = false;
		v->v.rr.distance.dist = std::numeric_limits<weight_t>::max();
		r.unbind_incoming(v->v.rr.node);

		if (v->v.circ.bias - v->v.circ.flow >= thresh) {
			TAILQ_INSERT_TAIL(&sources, v, v.rr.imbalance);
			source_cnt++;
			v->v.rr.instack = true;
			v->v.rr.distance.dist = 0;
			v->v.rr.distance.depth = 0;
			v->v.rr.distance.bias = 0;

			//printf("root: %d\n", v->id);

			RB_INSERT(vertextree, &stack, v);
		}
		else if (v->v.circ.bias - v->v.circ.flow <= -thresh)
			TAILQ_INSERT_TAIL(&sinks, v, v.rr.imbalance);
	}


	if (TAILQ_EMPTY(&sinks)) return;

	vertexhead workset, affected;

	while (true) {
		processstack(r, &stack);
		vertex_t *s = TAILQ_FIRST(&sinks);
		if (augmentpath(r, s, delta, thresh, &workset))
			source_cnt--;
		if (s->v.circ.bias - s->v.circ.flow > -thresh)
			TAILQ_REMOVE(&sinks, s, v.rr.imbalance);
		//printf("%d %d\n", source_cnt, s->id);
		if (TAILQ_EMPTY(&sinks) || source_cnt == 0) break;

		if (source_cnt % 10 == 0)
			fprintf(stderr, "imbalanced nodes: %d  \r", source_cnt);

		processworkset(r, &workset, &affected);
		processaffected(r, &affected, &stack);
	}

	//printf("\ndone %f\n", thresh);

}

void
canonize(network & g, residual & r, vertex_t *canon)
{
	vertextree stack;
	RB_INIT(&stack);

	for (vertex_t *v = g.get_first(); v; v = v->next()) {
		v->v.rr.instack = false;
		v->v.rr.inworkset = false;
		v->v.rr.isaffected = false;
		v->v.rr.completed = false;
		v->v.rr.distance.dist = std::numeric_limits<weight_t>::max();
		r.unbind_incoming(v->v.rr.node);
	}

	canon->v.rr.instack = true;
	canon->v.rr.distance.dist = 0;
	canon->v.rr.distance.depth = 0;
	canon->v.rr.distance.bias = 0;
	RB_INSERT(vertextree, &stack, canon);

	processstack(r, &stack);
	
	weight_t shift = canon->v.circ.dual;
	for (vertex_t *v = g.get_first(); v; v = v->next())
		v->v.circ.dual -= shift;
}


void
canonize(network & g, residual & r)
{
	vertextree stack;
	RB_INIT(&stack);

	for (vertex_t *v = g.get_first(); v; v = v->next()) {
		assert(v->v.circ.dual >= 0);
		v->v.rr.instack = true;
		v->v.rr.inworkset = false;
		v->v.rr.isaffected = false;
		v->v.rr.completed = false;
		v->v.rr.distance.dist = v->v.circ.dual;
		v->v.rr.distance.depth = 0;
		v->v.rr.distance.bias = 0;
		RB_INSERT(vertextree, &stack, v);

		r.unbind_incoming(v->v.rr.node);
	}

	processstack(r, &stack);
}
