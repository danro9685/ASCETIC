#include "circulation.h"
#include <stdio.h>

network *
read(FILE *f, bool weighted, bool self, uint32_t & cnt, uint32_t & ecnt, uint32_t limit)
{
	uint32_t a, b;

	cnt = 0;
	ecnt = 0;

	uintmap lm;


	while (fscanf(f, "%d%d", &a, &b) == 2) {
		if (weighted) fscanf(f, "%*f");
		if (!self && a == b) continue;
		if (lm.count(a) == 0)
			lm[a] = cnt++;
		if (lm.count(b) == 0)
			lm[b] = cnt++;
		ecnt++;
	}

	printf("%d vertices, %d edges\n", cnt, ecnt);

	rewind(f);

	network *g = new network(2*(cnt + ecnt + 2), 2*ecnt + 2*cnt + 1);
	for (uint32_t i = 0; i < cnt + ecnt + 2; i++) {
		vertex_t *v = g->addnode();
		v->v.init(v);
	}

	uint32_t ind = 0;
	while (fscanf(f, "%d%d", &a, &b) == 2) {
		double w = 1;
		if (weighted) fscanf(f, "%lf", &w);
		if (!self && a == b) continue;
		uint32_t x = lm[a];
		uint32_t y = lm[b];

		vertex_t *from = g->get(x);
		vertex_t *to = g->get(y);
		from->v.label = a;
		to->v.label = b;
		vertex_t *middle = g->get(ind + cnt);
		arc_t *e1 = g->addedge();
		arc_t *e2 = g->addedge();

		g->bindedge(e1, from, middle);
		g->bindedge(e2, to, middle);
		e2->v.circ.cost = 1;
		middle->v.circ.bias = -w;
		to->v.circ.bias += w;

		//printf("%d %d %d\n", from->id, to->id, middle->id);

		ind++;
	}

	vertex_t *src = g->get(cnt + ecnt);
	vertex_t *sink = g->get(cnt + ecnt + 1);
	for (uint32_t i = 0; i < cnt; i++) {
        arc_t *e1 = g->addedge();
		arc_t *e2 = g->addedge();

		vertex_t *v = g->get(i);

		g->bindedge(e1, src, v);
		g->bindedge(e2, v, sink);
	}

	arc_t *loop = g->addedge();
	g->bindedge(loop, sink, src);
	if (limit == 0)
		loop->v.circ.cost = cnt - 1;
	else
		loop->v.circ.cost = limit - 1;

	/*
	g = network_t(edges.begin(), edges.end(), ei.begin(), cnt);

	viterator vi = boost::vertices(g).first;

	for (intmap::iterator it = lm.begin(); it != lm.end(); ++it)
		g[vi[it->second]].label = it->first;
	
	for (uint32_t i = 0; i < cnt; i++) {
		g[vi[i]].id = i;
		g[vi[i]].rank = 0;
		g[vi[i]].ccomp = 0;
	}
	
	//boost::print_graph(g);
	*/

	return g;
}

compgraph *
read(FILE *f, bool weighted, bool self)
{
	uint32_t a, b;

	uint32_t cnt = 0;
	uint32_t ecnt = 0;

	uintmap lm;


	while (fscanf(f, "%d%d", &a, &b) == 2) {
		if (weighted) fscanf(f, "%*f");
		if (!self && a == b) continue;
		if (lm.count(a) == 0)
			lm[a] = cnt++;
		if (lm.count(b) == 0)
			lm[b] = cnt++;
		ecnt++;
	}

	printf("%d vertices, %d edges\n", cnt, ecnt);

	rewind(f);

	compgraph *g = new compgraph(cnt, ecnt); 
	for (uint32_t i = 0; i < cnt; i++) {
		g->addnode();
	}

	uint32_t ind = 0;
	while (fscanf(f, "%d%d", &a, &b) == 2) {
		double w = 1;
		if (weighted) fscanf(f, "%lf", &w);
		if (!self && a == b) continue;
		uint32_t x = lm[a];
		uint32_t y = lm[b];

		cvertex_t *from = g->get(x);
		cvertex_t *to = g->get(y);
		from->v.label = a;
		to->v.label = b;
		carc_t *e = g->addedge();

		g->bindedge(e, from, to);
		e->v.cost = w;

		ind++;
	}

	return g;
}
