#include "splittree.h"

void update(vertex_t *v, splitvertex *a, weight_t olddiff);

splitvertex *
init(vertexhead *h, svlist & tree)
{
	splitvertex *a = new splitvertex;

	a->negative = new vertexhead;
	a->negative_zero = new vertexhead;
	a->positive = new vertexhead;
	a->positive_zero = new vertexhead;
	TAILQ_INIT(a->positive_zero);
	TAILQ_INIT(a->positive);
	TAILQ_INIT(a->negative_zero);
	TAILQ_INIT(a->negative);

	vertex_t *v, *vnext;
	for (v = TAILQ_FIRST(h); v; v = vnext) {
		vnext = TAILQ_NEXT(v, v.group);
		TAILQ_INSERT_TAIL(a->negative, v, v.group);
		for (arc_t *e = v->first_out(); e; e = e->next_out()) {
			vertex_t *to = e->to;
			v->v.flux -= e->v.w;
			to->v.flux += e->v.w;
			v->v.deg++;
			to->v.deg++;
		}
	}

	// moves vertices to correct group
	for (v = TAILQ_FIRST(a->negative); v; v = vnext) {
		vnext = TAILQ_NEXT(v, v.group);
		update(v, a, 0);
	}

	tree.push_back(a);

	return a;
}

void
free(svlist & tree)
{
	for (svlist::iterator it = tree.begin(); it != tree.end(); ++it) {
		splitvertex *a = *it;
		delete a->negative;
		delete a->positive;
		delete a->negative_zero;
		delete a->positive_zero;
		delete a;
	}
}



bool
leftsmaller(vertexhead *s1, vertexhead *s2)
{
	uint32_t c1 = 0, c2 = 0;
	vertex_t *v1 = TAILQ_FIRST(s1);
	vertex_t *v2 = TAILQ_FIRST(s2);

	while ( !((v1 == NULL && c1 <= c2) || (v2 == NULL && c1 > c2))) {
		if (c1 <= c2) {
			c1 += v1->v.deg;
			v1 = TAILQ_NEXT(v1, v.group);
		}
		else {
			c2 += v2->v.deg;
			v2 = TAILQ_NEXT(v2, v.group);
		}
	}

	return v1 == NULL && c1 <= c2;
}


weight_t
gainleft(splitvertex *a)
{
	weight_t g = a->back + a->outback + a->diffneg;

	vertex_t *v;
	TAILQ_FOREACH(v, a->negative, v.group) {
		//printf("foo %f\n", g);
		g += v->v.diff();
	}
	//printf("foo %f\n", g);
	return g;	
}


weight_t
gainright(splitvertex *a)
{
	weight_t g = a->back + a->inback - a->diffpos;

	vertex_t *v;
	TAILQ_FOREACH(v, a->positive, v.group) {
		//printf("neg\n");
		g -= v->v.diff();
	}
	return g;	
}


void
split(splitvertex *a, network & g, svlist & tree)
{
	/*printf("LEAF: %p\n", a);
	vertex_t *v;
	TAILQ_FOREACH(v, a->negative, v.group) {
		printf("%d %f\n", v->v.label, v->v.diff());
	}
	TAILQ_FOREACH(v, a->positive, v.group) {
		printf("%d %f\n", v->v.label, v->v.diff());
	}
	TAILQ_FOREACH(v, a->negative_zero, v.group) {
		printf("%d %f\n", v->v.label, v->v.diff());
	}
	TAILQ_FOREACH(v, a->positive_zero, v.group) {
		printf("%d %f\n", v->v.label, v->v.diff());
	}
	*/

	if (leftsmaller(a->negative, a->positive)) {
		weight_t gain = gainleft(a);
		//printf("split left %f\n", gain);
		if (gain >= 0) return;
		assert(!TAILQ_EMPTY(a->negative));
		assert(!TAILQ_EMPTY(a->positive));
		a->gain = gain;
		splitleft(a, g, tree);
	}
	else {
		weight_t gain= gainright(a);
		//printf("split right %f\n", gain);
		if (gain >= 0) return;
		a->gain = gain;
		splitright(a, g, tree);
	}
}




void
update(vertex_t *v, splitvertex *a, weight_t olddiff)
{
	//printf("update: %d %p %f\n", v->v.label, a, olddiff);

	weight_t d = v->v.diff();

	if (!((d <= 0 && olddiff > 0) || (d > 0 && olddiff <= 0) || v->v.deg == 0)) 
		return;

	if (olddiff > 0)
		TAILQ_REMOVE(a->positive, v, v.group);
	if (olddiff <= 0)
		TAILQ_REMOVE(a->negative, v, v.group);

	if (v->v.deg == 0) {
		if (d <= 0) {
			TAILQ_INSERT_TAIL(a->negative_zero, v, v.group);
			a->diffneg += v->v.diff();
		}
		else {
			TAILQ_INSERT_TAIL(a->positive_zero, v, v.group);
			a->diffpos += v->v.diff();
		}
	}
	else {
		if (d <= 0)
			TAILQ_INSERT_TAIL(a->negative, v, v.group);
		else
			TAILQ_INSERT_TAIL(a->positive, v, v.group);
	}
}

void
splitleft(splitvertex *a, network & g, svlist & tree)
{
	splitvertex *left = new splitvertex;
	splitvertex *right = new splitvertex;

	a->left = left;
	a->right = right;

	//printf("split: %p %p\n", left, right);


	left->negative = a->negative;
	left->negative_zero = a->negative_zero;
	left->positive = new vertexhead;
	TAILQ_INIT(left->positive);
	left->positive_zero = new vertexhead;
	TAILQ_INIT(left->positive_zero);

	right->positive = a->positive;
	right->positive_zero = a->positive_zero;
	right->negative = new vertexhead;
	TAILQ_INIT(right->negative);
	right->negative_zero = new vertexhead;
	TAILQ_INIT(right->negative_zero);

	left->back = a->back + a->outback;
	left->inback = 0;
	left->outback = 0;
	left->diffneg = a->diffneg;
	left->diffpos = 0;


	right->back = a->back;
	right->inback = a->inback;
	right->outback = a->outback;
	right->diffneg = 0;
	right->diffpos = a->diffpos;

	a->negative = 0;
	a->positive = 0;
	a->negative_zero = 0;
	a->positive_zero = 0;


	vertex_t *v, *vnext;
	TAILQ_FOREACH(v, left->negative, v.group) {
		arc_t *e;
		for (e = v->first_in(); e; e = e->next_in()) {
			vertex_t *w = e->from;
			if (w->v.diff() > 0) e->v.cross = true;
		}
		for (e = v->first_out(); e; e = e->next_out()) {
			vertex_t *w = e->to;
			if (w->v.diff() > 0) e->v.cross = true; 
		}
	}

	for (v = TAILQ_FIRST(left->negative); v; v = vnext) {
		vnext = TAILQ_NEXT(v, v.group);

		left->back -= v->v.outback;
		right->back += v->v.inback;

		left->outback += v->v.outback;
		right->outback -= v->v.outback;

		left->inback += v->v.inback;
		right->inback -= v->v.inback;

		arc_t *e, *enext;

		for (e = v->first_in(); e; e = enext) {
			enext = e->next_in();
			if (!e->v.cross) continue; // not a cross edge

			//printf("deleting in\n");

			vertex_t *w = e->from;

			weight_t vdiff = v->v.diff();
			weight_t wdiff = w->v.diff();
		
			g.unbindedge(e);
			v->v.deg--;
			v->v.flux -= e->v.w;
			v->v.inback += e->v.w;
			left->inback += e->v.w;

			w->v.deg--;
			w->v.flux += e->v.w;
			w->v.outback += e->v.w;
			right->outback += e->v.w;

			update(v, left, vdiff);
			update(w, right, wdiff);
		}

		for (e = v->first_out(); e; e = enext) {
			enext = e->next_out();
			if (!e->v.cross) continue; // not a cross edge
			vertex_t *w = e->to;

			//printf("deleting out\n");


			weight_t vdiff = v->v.diff();
			weight_t wdiff = w->v.diff();
		
			g.unbindedge(e);
			v->v.deg--;
			v->v.flux += e->v.w;

			w->v.deg--;
			w->v.flux -= e->v.w;

			update(v, left, vdiff);
			update(w, right, wdiff);
		}
	}
	tree.push_back(left);
	tree.push_back(right);
	split(left, g, tree);
	split(right, g, tree);
}


void
splitright(splitvertex *a, network & g, svlist & tree)
{
	splitvertex *left = new splitvertex;
	splitvertex *right = new splitvertex;

	a->left = left;
	a->right = right;

	//printf("split: %p %p\n", left, right);

	left->negative = a->negative;
	left->negative_zero = a->negative_zero;
	left->positive = new vertexhead;
	TAILQ_INIT(left->positive);
	left->positive_zero = new vertexhead;
	TAILQ_INIT(left->positive_zero);

	right->positive = a->positive;
	right->positive_zero = a->positive_zero;
	right->negative = new vertexhead;
	TAILQ_INIT(right->negative);
	right->negative_zero = new vertexhead;
	TAILQ_INIT(right->negative_zero);

	left->back = a->back;
	left->inback = a->inback;
	left->outback = a->outback;
	left->diffneg = a->diffneg;
	left->diffpos = 0;


	right->back = a->back + a->inback;
	right->inback = 0;
	right->outback = 0;
	right->diffneg = 0;
	right->diffpos = a->diffpos;

	a->negative = 0;
	a->positive = 0;
	a->negative_zero = 0;
	a->positive_zero = 0;

	vertex_t *v, *vnext;

	//mark cross edges
	TAILQ_FOREACH(v, right->positive, v.group) {
		arc_t *e;
		for (e = v->first_in(); e; e = e->next_in()) {
			vertex_t *w = e->from;
			if (w->v.diff() <= 0) e->v.cross = true;
		}
		for (e = v->first_out(); e; e = e->next_out()) {
			vertex_t *w = e->to;
			if (w->v.diff() <= 0) e->v.cross = true; 
		}
	}


	for (v = TAILQ_FIRST(right->positive); v; v = vnext) {
		vnext = TAILQ_NEXT(v, v.group);

		left->back += v->v.outback;
		right->back -= v->v.inback;

		left->outback -= v->v.outback;
		right->outback += v->v.outback;

		left->inback -= v->v.inback;
		right->inback += v->v.inback;

		arc_t *e, *enext;
		for (e = v->first_out(); e; e = enext) {
			enext = e->next_out();
			if (!e->v.cross) continue; // not a cross edge

			vertex_t *w = e->to;

			weight_t vdiff = v->v.diff();
			weight_t wdiff = w->v.diff();

			g.unbindedge(e);
			v->v.deg--;
			v->v.flux += e->v.w;
			v->v.outback += e->v.w;
			right->outback += e->v.w;

			w->v.deg--;
			w->v.flux -= e->v.w;
			w->v.inback += e->v.w;
			left->inback += e->v.w;


			update(v, right, vdiff);
			update(w, left, wdiff);
		}

		for (e = v->first_in(); e; e = enext) {
			enext = e->next_in();
			if (!e->v.cross) continue; // not a cross edge

			vertex_t *w = e->from;

			weight_t vdiff = v->v.diff();
			weight_t wdiff = w->v.diff();
		
			g.unbindedge(e);
			v->v.deg--;
			v->v.flux -= e->v.w;

			w->v.deg--;
			w->v.flux += e->v.w;

			update(v, right, vdiff);
			update(w, left, wdiff);
		}
	}

	tree.push_back(left);
	tree.push_back(right);
	split(left, g, tree);
	split(right, g, tree);
}
