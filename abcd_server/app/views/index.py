from flask import Blueprint, render_template, flash, redirect, url_for

bp = Blueprint('index', __name__)

index_title = {
    'title': 'Welcome',
    'subtitle': 'ABCD database'
}

index_content = """
Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the
standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to
make a type specimen book. It has survived not only five centuries, but also the leap into electronic
typesetting, remaining essentially unchanged. It was popularised in the 1960s with the release of Letraset
sheets containing Lorem Ipsum passages, and more recently with desktop publishing software like Aldus
PageMaker including versions of Lorem Ipsum.
"""


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route('/')
def index():
    return render_template("index.html",
                           title=index_title,
                           content=index_content)

# from flask_bootstrap import __version__ as FLASK_BOOTSTRAP_VERSION
# from markupsafe import escape
# from forms import SignupForm

# from app.nav import nav
# from database import generate_plot, Item, ItemTable
# from admin import requires_auth

# @frontend.app_template_filter('reverse')
# def include_file(name):
#     return name+'haha'

# # Shows a long signup form, demonstrating form rendering.
# @frontend.route('/example-form/', methods=('GET', 'POST'))
# def example_form():
#     form = SignupForm()
#
#     if form.validate_on_submit():
#         # We don't have anything fancy in our application, so we are just
#         # flashing a message when a user completes the form successfully.
#         #
#         # Note that the default flashed messages rendering allows HTML, so
#         # we need to escape things if we input user values:
#         flash('Hello, {}. You have successfully signed up'
#               .format(escape(form.name.data)))
#
#         # In a real application, you may wish to avoid this tedious redirect.
#         return redirect(url_for('.index'))
#
#     return render_template('signup.html', form=form)
