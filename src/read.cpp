#include "circulation.h"
#include <stdio.h>
#include <Rcpp.h>
using namespace Rcpp;

network *
read(FILE *f, bool weighted, bool self, uint32_t & cnt, uint32_t & ecnt, uint32_t limit)
{
	uint32_t a, b;

	cnt = 0;
	ecnt = 0;

	uintmap lm;


	while (fscanf(f, "%d%d", &a, &b) == 2) {
		if (weighted) {
			int resscanf = fscanf(f, "%*f");
		}
		if (!self && a == b) {
			continue;
		}
		if (lm.count(a) == 0) {
			lm[a] = cnt++;
		}
		if (lm.count(b) == 0) {
			lm[b] = cnt++;
		}
		ecnt++;
	}

	//printf("%d vertices, %d edges\n", cnt, ecnt);

	rewind(f);

	network *g = new network(2*(cnt + ecnt + 2), 2*ecnt + 2*cnt + 1);
	for (uint32_t i = 0; i < cnt + ecnt + 2; i++) {
		vertex_t *v = g->addnode();
		v->v.init(v);
	}

	uint32_t ind = 0;
	while (fscanf(f, "%d%d", &a, &b) == 2) {
		double w = 1;
		if (weighted) {
			int resscanf = fscanf(f, "%lf", &w);
		}
		
		if (!self && a == b) {
			continue;
		}
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
read(FILE *f, IntegerMatrix m, bool weighted, bool self)
{
	uint32_t a, b;
  uint32_t a_alt, b_alt;

	uint32_t cnt = 0;
	uint32_t ecnt = 0;
	
	uint32_t cnt_alt = 0;
	uint32_t ecnt_alt = 0;

	uintmap lm;
	uintmap lm_alt;

  // to be rewritten in a for loop, check what to do with "weighted"
  // resscanf is an integer, assuming third possible column type is integer
  
  for(int i = 0; i < m.nrow(); i++) {
    a = m(i,0);
    b = m(i,1);
    printf("a: %d - b: %d \n", a, b);
    printf("wighted: %s \n", (weighted ? "true" : "false"));
    if (weighted) {
      int resscanf = m(i,2);
      printf("resscanf: %u\n", resscanf);
    }
    
    if (a != b || self) {
      if (lm.count(a) == 0) {
        lm[a] = cnt++;
      }
      if (lm.count(b) == 0) {
        lm[b] = cnt++;
      }
      ecnt++;
    }
    
  }
  
  
	while (fscanf(f, "%d%d", &a_alt, &b_alt) == 2) {
	  
	  printf("a: %d - b: %d \n", a_alt, b_alt);
	  
	  printf("wighted: %s \n", (weighted ? "true" : "false"));
	  
		if (weighted) {
			int resscanf = fscanf(f, "%*f");
		  printf("resscanf: %u\n", resscanf);
		}
		
		if (!self && a_alt == b_alt) {
			continue;
		}
		if (lm_alt.count(a_alt) == 0) {
			lm_alt[a_alt] = cnt_alt++;
		}
		if (lm_alt.count(b_alt) == 0) {
			lm_alt[b_alt] = cnt_alt++;
		}
		ecnt_alt++;
	}

	printf("%d vertices, %d edges NORMAL\n", cnt, ecnt);
  printf("%d vertices, %d edges ALT\n", cnt_alt, ecnt_alt);

	rewind(f);

	compgraph *g = new compgraph(cnt, ecnt); 
	for (uint32_t i = 0; i < cnt; i++) {
		g->addnode();
	}
	
	compgraph *g_alt = new compgraph(cnt_alt, ecnt_alt); 
	for (uint32_t i = 0; i < cnt_alt; i++) {
	  g_alt->addnode();
	}

	
	// ok now this other loop
	uint32_t ind = 0;
	uint32_t ind_alt = 0;
	
	for(int i = 0; i < m.nrow(); i++) {
	  a = m(i,0);
	  b = m(i,1);
	  printf("a: %d - b: %d \n", a, b);
	  printf("wighted: %s \n", (weighted ? "true" : "false"));
	  double w = 1;
	  
	  if (weighted) {
	    int resscanf = m(i,2);
	    printf("resscanf: %u\n", resscanf);
	  }
	  
	  if (a != b || self) {
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
	  
	}
	
	while (fscanf(f, "%d%d", &a_alt, &b_alt) == 2) {
	  printf("a: %d - b: %d \n", a_alt, b_alt);
		double w = 1;
		if (weighted) {
			int resscanf = fscanf(f, "%lf", &w);
		  printf("resscanf: %u\n", resscanf);
		}
		
		if (!self && a_alt == b_alt) {
			continue;
		}
		uint32_t x = lm_alt[a_alt];
		uint32_t y = lm_alt[b_alt];

		cvertex_t *from = g_alt->get(x);
		cvertex_t *to = g_alt->get(y);
		from->v.label = a_alt;
		to->v.label = b_alt;
		carc_t *e = g_alt->addedge();

		g_alt->bindedge(e, from, to);
		e->v.cost = w;

		ind_alt++;
	}
	
	printf("ind %u\n", ind);
	printf("ind_alt %u\n", ind_alt);

	return g;
}
