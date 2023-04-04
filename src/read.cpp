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

	return g;
}

compgraph * read(IntegerMatrix m)
{
	uint32_t a, b;

	uint32_t cnt = 0;
	uint32_t ecnt = 0;

	uintmap lm;

  for(int i = 0; i < m.nrow(); i++) {
    a = m(i,0);
    b = m(i,1);
    
    if (a != b) {
      if (lm.count(a) == 0) {
        lm[a] = cnt++;
      }
      if (lm.count(b) == 0) {
        lm[b] = cnt++;
      }
      ecnt++;
    }
    
  }
  
	compgraph *g = new compgraph(cnt, ecnt); 
	for (uint32_t i = 0; i < cnt; i++) {
		g->addnode();
	}
	
	uint32_t ind = 0;
	
	for(int i = 0; i < m.nrow(); i++) {
	  a = m(i,0);
	  b = m(i,1);
	  double w = 1;
	  
	  if (a != b) {
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

	return g;
}
