#include "circulation.h"


TAILQ_HEAD(carchead, carc_t);

uint32_t
ccomp(compgraph & g, cvertexhead *roots)
{
	carchead edges;
	cvertexhead stack;
	uint32_t index = 1;
	uint32_t compid = 0;

	TAILQ_INIT(&edges);
	TAILQ_INIT(&stack);
	TAILQ_INIT(roots);

	for (cvertex_t *r = g.get_first(); r; r = r->next()) {
		if (r->v.index > 0) continue;

		cvertex_t *v = r;
		carc_t *e = 0;


		while (true) {
			if (e == 0 && v->v.index == 0) {
				v->v.index = index;
				v->v.lowlink = index;
				index++;
				// push to stack
				TAILQ_INSERT_HEAD(&stack, v, v.stack);
				v->v.instack = true;
				e = v->first_out();
			}
			else if (e == 0 && v->v.index > 0) {
				if (v->v.lowlink == v->v.index) {
					// scc found
					cvertex_t *w;
					uint32_t cid = 0;
					TAILQ_INIT(&v->v.comps);
					TAILQ_INSERT_TAIL(roots, v, v.root);
					do {
						w = TAILQ_FIRST(&stack);
						TAILQ_REMOVE(&stack, w, v.stack);
						w->v.compid = compid;
						w->v.compind = cid++;
						w->v.instack = false;
						//printf("%d %d\n", w->id, compid);
						TAILQ_INSERT_HEAD(&v->v.comps, w, v.comp);
					} while (v != w);
					compid++;
				}
				e = TAILQ_FIRST(&edges);
				if (e == 0) break;
				TAILQ_REMOVE(&edges, e, v.stack);
				v = e->from;
				v->v.lowlink = std::min(v->v.lowlink, e->to->v.lowlink);
				e = e->next_out();
			}
			else {
				cvertex_t *w = e->to;
				if (w->v.index == 0) {
					TAILQ_INSERT_HEAD(&edges, e, v.stack);
					v = w;
					e = 0;
				}
				else {
					if (w->v.instack)
						v->v.lowlink = std::min(v->v.lowlink, w->v.lowlink);
					e = e->next_out();
				}
			}
		}

	}
	return compid;
}

network *
create_network(compgraph & g, cvertex_t *root, uint32_t & ncnt, uint32_t & ecnt, uint32_t limit)
{
	uint32_t compid = root->v.compid;
	ncnt = 0;
	ecnt = 0;

	// compute node and edge size
	cvertex_t *v;
	TAILQ_FOREACH(v, &root->v.comps, v.comp) {
		if (v->v.compid == compid) {
			ncnt++;
			for (carc_t *e = v->first_out(); e; e = e->next_out()) {
				if (e->to->v.compid == compid) ecnt++;
			}
		}
	}

	network *n = new network(2*(ncnt + ecnt + 2), 2*ecnt + 2*ncnt + 1);

	for (uint32_t i = 0; i < ncnt + ecnt + 2; i++) {
		vertex_t *v = n->addnode();
		v->v.init(v);
	}
	
	uint32_t ind = 0;
	TAILQ_FOREACH(v, &root->v.comps, v.comp) {
		if (v->v.compid == compid) {
			for (carc_t *e = v->first_out(); e; e = e->next_out()) {
				if (e->to->v.compid == compid) {
					vertex_t *from = n->get(v->v.compind);
					vertex_t *to = n->get(e->to->v.compind);
					from->v.label = 0;
					to->v.label = 0;
					vertex_t *middle = n->get(ind + ncnt);
					arc_t *e1 = n->addedge();
					arc_t *e2 = n->addedge();

					n->bindedge(e1, from, middle);
					n->bindedge(e2, to, middle);
					e2->v.circ.cost = 1;
					middle->v.circ.bias = -e->v.cost;
					to->v.circ.bias += e->v.cost;

					//printf("%f\n", e->v.cost);

					ind++;
				}
			}
		}
	}


	vertex_t *src = n->get(ncnt + ecnt);
	vertex_t *sink = n->get(ncnt + ecnt + 1);
	for (uint32_t i = 0; i < ncnt; i++) {
		arc_t *e1 = n->addedge();
		arc_t *e2 = n->addedge();

		vertex_t *v = n->get(i);

		n->bindedge(e1, src, v);
		n->bindedge(e2, v, sink);
	}

	arc_t *loop = n->addedge();
	n->bindedge(loop, sink, src);
	if (limit == 0)
		loop->v.circ.cost = ncnt - 1;
	else
		loop->v.circ.cost = limit - 1;

	return n;
}

// joins the connected components into one big graph so that we can get a canonical solution.
network *
join(networkvector & comps, compgraph & cg)
{
	uint32_t ncnt = 0;
	uint32_t ecnt = 0;

	uintvector shift(comps.size());
	for (int32_t i = comps.size() - 1; i >= 0; i--) {
		network *g = comps[i];
		shift[i] = ncnt;
		for (vertex_t *v = g->get_first(); v; v = v->next()) {
			v->v.label = ncnt++;
			for (arc_t *e = v->first_out(); e; e = e->next_out()) ecnt++;
		}
	}



	network *j = new network(ncnt, ecnt + cg.edgebudget()); 

	for (int32_t i = comps.size() - 1; i >= 0; i--) {
		network *g = comps[i];
		// copy nodes
		for (vertex_t *v = g->get_first(); v; v = v->next()) {
			vertex_t *w = j->addnode();
			w->v.init(w);
			w->v.circ.dual = v->v.circ.dual + shift[i];
		}
	}

	for (uint32_t i = 0; i < comps.size(); i++) {
		network *g = comps[i];
		// copy edges
		for (vertex_t *v = g->get_first(); v; v = v->next()) {
			for (arc_t *e = v->first_out(); e; e = e->next_out()) {
				vertex_t *w = e->to;
				arc_t *f = j->addedge();
				j->bindedge(f, j->get(v->v.label), j->get(w->v.label));
				f->v.circ.flow = e->v.circ.flow;
				f->v.circ.cost = e->v.circ.cost;
			}
		}
		
	}

	// add cross-component edges
	for (cvertex_t *v = cg.get_first(); v; v = v->next()) {
		for (carc_t *e = v->first_out(); e; e = e->next_out()) {
			cvertex_t *w = e->to;
			if (v->v.compid == w->v.compid) continue;

			vertex_t *x = comps[v->v.compid]->get(v->v.compind);
			vertex_t *y = comps[w->v.compid]->get(w->v.compind);
			vertex_t *mx = x->v.contract.master;
			vertex_t *my = y->v.contract.master;

			vertex_t *a = j->get(mx->v.label);
			vertex_t *b = j->get(my->v.label);
			arc_t *f = j->addedge();
			j->bindedge(f, a, b);

			f->v.circ.cost = -1 - (x->v.circ.dual - mx->v.circ.dual) + (y->v.circ.dual - my->v.circ.dual);
			//printf("%d %d %d %d %d %d %d\n", f->v.circ.cost, a->id, b->id, v->v.compid, w->v.compid, v->id, w->id);
		}
	}

	return j;
}

void
migrate_dual(compgraph & cg, networkvector & comps, network & j)
{
	for (cvertex_t *v = cg.get_first(); v; v = v->next()) {
		vertex_t *x = comps[v->v.compid]->get(v->v.compind);
		vertex_t *mx = x->v.contract.master;
		vertex_t *a = j.get(mx->v.label);

		v->v.dual = a->v.circ.dual + (x->v.circ.dual - mx->v.circ.dual);

	}
}
