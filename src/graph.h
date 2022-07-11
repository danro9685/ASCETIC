#ifndef GRAPH_H
#define GRAPH_H

#include "queue.h"
#include "defines.h"
#include <assert.h>

template<typename NT, typename ET> struct edge_t;

template<typename NT, typename ET>
struct node_t {
			typedef node_t<NT, ET> ntype;
			typedef edge_t<NT, ET> etype;

			NT v;
			uint32_t id;

			TAILQ_HEAD(edgehead, etype) in, out;
			TAILQ_ENTRY(ntype) active;

			etype *first_out() {return TAILQ_FIRST(&out);}
			etype *first_in() {return TAILQ_FIRST(&in);}

			TAILQ_HEAD(nodehead, ntype);


			ntype *next() {return TAILQ_NEXT(this, active);}
			ntype *prev() {return TAILQ_PREV(this, nodehead, active);}

		};

template<typename NT, typename ET>
		struct edge_t {
			typedef node_t<NT, ET> ntype;
			typedef edge_t<NT, ET> etype;

			ET v;
			ntype *from, *to;
			TAILQ_ENTRY(etype) in, out;
			
			etype *next_out() {return TAILQ_NEXT(this, out);}
			etype *next_in() {return TAILQ_NEXT(this, in);}
		};




template <typename NT, typename ET>
class graph {
	public:
		typedef node_t<NT, ET> node;
		typedef edge_t<NT, ET> edge;

		TAILQ_HEAD(nodehead, node);


		graph(uint32_t n, uint32_t m) :
			m_nodes(n), m_edges(m),
			m_nodefree(0), m_edgefree(0), m_activenodes(0)
		{
			TAILQ_INIT(&m_active);
			for (uint32_t i = 0; i < n; i++) m_nodes[i].id = i;
		}

		node *
		addnode()
		{
			node *n = &m_nodes[m_nodefree++];
			TAILQ_INIT(&n->in);
			TAILQ_INIT(&n->out);
			TAILQ_INSERT_TAIL(&m_active, n, active);
			m_activenodes++;
			return n;
		}

		void
		deactivatenode(node *n)
		{
			TAILQ_REMOVE(&m_active, n, active);
			m_activenodes--;
		}

		edge *
		addedge()
		{
			edge *e = &m_edges[m_edgefree++];
			e->from = e->to = 0;
			return e;
		}

		void 
		unbindedge(edge *e)
		{
			if (e->from) {
				TAILQ_REMOVE(&e->from->out, e, out);
				TAILQ_REMOVE(&e->to->in, e, in);
				e->from = e->to = 0;
			}
		}

		void
		bindedge(edge *e, node *from, node *to)
		{
			assert(from != to);
			unbindedge(e);
			e->from = from;
			e->to = to;

			TAILQ_INSERT_TAIL(&from->out, e, out);
			TAILQ_INSERT_TAIL(&to->in, e, in);
		}

		void
		unbind_incoming(node *to)
		{
			edge *e;
			while ((e = to->first_in()) != 0)
				unbindedge(e);
		}

		typedef std::vector<node *> nodevector;
		typedef std::vector<edge *> edgevector;

		node *get_first() {return TAILQ_FIRST(&m_active);}
		node *get_last() {return TAILQ_LAST(&m_active, nodehead);}
		node *get(uint32_t i) {return &m_nodes[i];}

		uint32_t nodebudget() const {return m_nodes.size();}
		uint32_t edgebudget() const {return m_edges.size();}
		uint32_t nodecnt() const {return m_activenodes;}

		void
		wipe()
		{
			TAILQ_INIT(&m_active);
			m_nodefree = 0;
			m_edgefree = 0;
			m_activenodes = 0;
		}

	protected:
		typedef std::vector<node> nodearray;
		typedef std::vector<edge> edgearray;

		nodearray m_nodes;
		edgearray m_edges;

		uint32_t m_nodefree, m_edgefree, m_activenodes;

		nodehead m_active;
};

#endif
