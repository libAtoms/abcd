from flask_nav import Nav
from flask_nav.elements import Navbar, View, Separator, Subgroup, Link
from hashlib import sha1

from dominate import tags
from flask_nav.renderers import Renderer

nav = Nav()


@nav.navigation()
def topnavbar():
    return Navbar(
        'ABCD',
        View('Database', 'index.database'),
        View('GraphQL', 'index.schema'),
        Subgroup(
            'Links',
            Link('Docs', 'https://fekad.github.io/abcd/'),
            Separator(),
            Link('Github', 'https://github.com/fekad/abcd'),
        ),
    )


class BootstrapRenderer(Renderer):
    def __init__(self, nav_id=None):
        self.id = nav_id

        self._in_dropdown = False

    def visit_Navbar(self, node):
        node_id = self.id or sha1(str(id(node)).encode()).hexdigest()

        root = tags.nav(_class='navbar navbar-expand-md navbar-dark bg-dark fixed-top')

        # title may also have a 'get_url()' method, in which case we render
        # a brand-link
        if node.title is not None:
            if hasattr(node.title, 'get_url'):
                root.add(tags.a(node.title.text, _class='navbar-brand',
                                href=node.title.get_url()))
            else:
                root.add(tags.span(node.title, _class='navbar-brand'))

        btn = root.add(tags.button())
        btn['type'] = 'button'
        btn['class'] = 'navbar-toggler'
        btn['data-toggle'] = 'collapse'
        btn['data-target'] = '#' + node_id
        btn['aria-controls'] = 'navbarCollapse'
        btn['aria-expanded'] = 'false'
        btn['aria-label'] = "Toggle navigation"

        btn.add(tags.span('', _class='navbar-toggler-icon'))

        bar = root.add(tags.div(
            _class='navbar-collapse collapse',
            id=node_id,
        ))

        bar_list = bar.add(tags.ul(_class='navbar-nav mr-auto'))

        for item in node.items:
            bar_list.add(self.visit(item))

        # search_form = root.add(tags.form(_class="form-inline mt-2 mt-md-0"))
        # search_input = search_form.add(tags.input(_class="form-control mr-sm-2"))
        # search_input['type'] = "text"
        # search_input['placeholder'] = "Search"
        # search_input['aria-label'] = "Search"
        #
        # search_btn = search_form.add(tags.button(_class="btn btn-outline-success my-2 my-sm-0"))
        # search_btn['type'] = "submit"
        # search_btn.add_raw_string('Search')

        return root

    def visit_Text(self, node):
        if self._in_dropdown:
            return tags.a(node.text, _class='dropdown-item disabled')

        return tags.li(tags.a(node.text, _class='nav-link disabled'), _class="nav-item")

    def visit_Link(self, node):
        if self._in_dropdown:
            return tags.a(node.text, href=node.get_url(), _class='dropdown-item')

        item = tags.li(_class='nav-item')
        item.add(tags.a(node.text, href=node.get_url(), _class='nav-link'))

        return item

    def visit_View(self, node):
        if self._in_dropdown:
            return tags.a(node.text, href=node.get_url(), _class='dropdown-item')

        item = tags.li(_class='nav-item')
        item.add(tags.a(node.text, href=node.get_url(), title=node.text, _class='nav-link'))
        if node.active:
            item['class'] = 'nav-item active'

        return item

    def visit_Subgroup(self, node):

        if self._in_dropdown:
            raise RuntimeError('Cannot render nested Subgroups')

        li = tags.li(_class='nav-item dropdown')
        if node.active:
            li['class'] = 'nav-item dropdown active'

        a = li.add(tags.a(node.title, href='#', _class='nav-link dropdown-toggle'))
        a['data-toggle'] = 'dropdown'
        a['aria-haspopup'] = 'true'
        a['aria-expanded'] = 'false'

        dropdown_div = li.add(tags.div(_class='dropdown-menu'))

        self._in_dropdown = True
        for item in node.items:
            dropdown_div.add(self.visit(item))
        self._in_dropdown = False

        return li

    def visit_Separator(self, node):
        if self._in_dropdown:
            return tags.div(_class='dropdown-divider')

        raise RuntimeError('Cannot render separator outside Subgroup.')
