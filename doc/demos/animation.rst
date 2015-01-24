Animations with graph-tool
==========================

The drawing capabilities of ``graph-tool`` (see :mod:`~graph_tool.draw`
module) can be harnessed to perform animations in a straightforward
manner. Here we show some examples which uses `GTK+
<http://www.gtk.org/>`_ to display animations in an
:class:`~graph_tool.draw.interactive_window`, as well as offscreen to a
file. The idea is to easily generate visualisations which can be used in
presentations, and embedded in websites.


SIRS epidemics
--------------

Here we implement a simple `SIRS epidemics
<http://en.wikipedia.org/wiki/Epidemic_model>`_ on a network, and we
construct an animation showing the time evolution. Nodes which are
susceptible (S) are shown in white, whereas infected (I) nodes are shown
in black. Recovered (R) nodes are removed from the layout, since they
cannot propagate the outbreak.

The script which performs the animation is called
:download:`animation_sirs.py <animation_sirs.py>` and is shown below.

.. literalinclude:: animation_sirs.py
   :linenos:


If called without arguments, the script will show the animation inside an
:class:`~graph_tool.draw.interactive_window`. If the parameter
``offscreen`` is passed, individual frames will be saved in the
``frames`` directory:

.. code-block:: bash

   $ ./animation_sirs.py offscreen

.. doctest::
   :hide:

   >>> import subprocess
   >>> subprocess.call(["demos/animation_sirs.py", "offscreen"])
   0

These frames can be combined and encoded into the appropriate
format. Here we use the `mencoder
<http://www.mplayerhq.hu/DOCS/HTML/en/mencoder.html>`_ tool from
`mplayer <http://www.mplayerhq.hu>`_ to combine all the frames into a
single file with YUY format, and then we encode this with the `WebM
format <http://www.webmproject.org>`_, using `vpxenc
<http://www.webmproject.org/docs/encoder-parameters/>`_, so that it can
be embedded in a website.

.. code-block:: bash

   $ mencoder mf://frames/sirs*.png -mf w=500:h=400:type=png -ovc raw -of rawvideo -vf format=i420 -nosound -o sirs.yuy
   $ vpxenc sirs.yuy -o sirs.webm -w 500 -h 400 --fps=25/1 --target-bitrate=1000 --good --threads=4

.. doctest::
   :hide:

   >>> import subprocess
   >>> subprocess.call("mencoder mf://frames/sirs*.png -mf w=500:h=400:type=png -ovc raw -of rawvideo -vf format=i420 -nosound -o demos/sirs.yuy".split())
   0
   >>> subprocess.call("vpxenc demos/sirs.yuy -o demos/sirs.webm -w 500 -h 400 --fps=25/1 --target-bitrate=1000 --good --threads=4".split())
   0


The resulting animation can be downloaded :download:`here <sirs.webm>`,
or played below if your browser supports WebM.

.. raw:: html

   <div style="text-align:center">
       <video id="sirs" src="../_downloads/sirs.webm" controls></video>
   </div>


This type of animation can be extended or customized in many ways, by
dynamically modifying the various drawing parameters and vertex/edge
properties. For instance, one might want to represent the susceptible
state as either |susceptible| or |susceptible-fear|, depending on
whether a neighbor is infected, and the infected state as |zombie|.
Properly modifying the script above would lead to the following
:download:`movie <zombie.webm>`:

.. doctest::
   :hide:

   >>> import subprocess
   >>> subprocess.call(["demos/animation_zombies.py", "offscreen"])
   0
   >>> import subprocess
   >>> subprocess.call("mencoder mf://frames/zombies*.png -mf w=500:h=400:type=png -ovc raw -of rawvideo -vf format=i420 -nosound -o demos/zombie.yuy".split())
   0
   >>> subprocess.call("vpxenc demos/zombie.yuy -o demos/zombie.webm -w 500 -h 400 --fps=10/1 --target-bitrate=1000 --good --threads=4".split())
   0

.. raw:: html

   <div style="text-align:center">
       <video id="sirs" src="../_downloads/zombie.webm" controls></video>
   </div>

The modified script can be downloaded :download:`here <animation_zombies.py>`.


.. |susceptible| image:: face-grin.png
   :height: 48
   :width: 48
.. |susceptible-fear| image:: face-surprise.png
   :height: 48
   :width: 48
.. |zombie| image:: zombie.png
   :height: 48
   :width: 48



Dynamic layout
--------------

The graph layout can also be updated during an animation. As an
illustration, here we consider a very simplistic model for spatial
segregation, where the edges of the graph are repeatedly and randomly
rewired, as long as the new edge has a shorter euclidean distance.

The script which performs the animation is called
:download:`animation_dancing.py <animation_dancing.py>` and is shown below.

.. literalinclude:: animation_dancing.py
   :linenos:


This example works like the SIRS example above, and if we pass the
``offscreen`` parameter, the frames will be dumped to disk, otherwise
the animation is displayed inside an :class:`~graph_tool.draw.interactive_window`.

.. code-block:: bash

   $ ./animation_dancing.py offscreen

.. doctest::
   :hide:

   >>> import subprocess
   >>> subprocess.call(["demos/animation_dancing.py", "offscreen"])
   0


Also like the previous example, we can encode the animation with the `WebM
format <http://www.webmproject.org>`_:

.. code-block:: bash

   $ mencoder mf://frames/dancing*.png -mf w=500:h=400:type=png -ovc raw -of rawvideo -vf format=i420 -nosound -o dancing.yuy
   $ vpxenc dancing.yuy -o dancing.webm -w 500 -h 400 --fps=100/1 --target-bitrate=5000 --good --threads=4


.. doctest::
   :hide:

   >>> import subprocess
   >>> subprocess.call("mencoder mf://frames/dancing*.png -mf w=500:h=400:type=png -ovc raw -of rawvideo -vf format=i420 -nosound -o demos/dancing.yuy".split())
   0
   >>> subprocess.call("vpxenc demos/dancing.yuy -o demos/dancing.webm -w 500 -h 400 --fps=100/1 --target-bitrate=2000 --good --threads=4".split())
   0


The resulting animation can be downloaded :download:`here
<dancing.webm>`, or played below if your browser supports WebM.

.. raw:: html

   <div style="text-align:center">
       <video id="sirs" src="../_downloads/dancing.webm" controls></video>
   </div>

Interactive visualizations
--------------------------

Here we show an example of interactive visualization where the BFS tree
of the currently selected vertex is highlighted with a different color.

The script which performs the visualization is called
:download:`interactive_bst.py <interactive_bst.py>` and is shown
below. When called, it will open an interactive window.

.. literalinclude:: interactive_bst.py
   :linenos:

.. code-block:: bash

   $ ./interactive_bst.py offscreen

.. doctest::
   :hide:

   >>> import subprocess
   >>> subprocess.call(["demos/interactive_bst.py", "offscreen"])
   0

The above script is interactive, i.e. it expects a reaction from the
user. But for the purpose of this demo, it also saves the frames to a
file, so we can encode the animation with the `WebM format
<http://www.webmproject.org>`_:

.. code-block:: bash

   $ mencoder mf://frames/bfs*.png -mf w=500:h=400:type=png -ovc raw -of rawvideo -vf format=i420 -nosound -o bfs.yuy
   $ vpxenc bfs.yuy -o bfs.webm -w 500 -h 400 --fps=5/1 --target-bitrate=5000 --good --threads=4


.. doctest::
   :hide:

   >>> import subprocess
   >>> subprocess.call("mencoder mf://frames/bfs*.png -mf w=500:h=400:type=png -ovc raw -of rawvideo -vf format=i420 -nosound -o demos/bfs.yuy".split())
   0
   >>> subprocess.call("vpxenc demos/bfs.yuy -o demos/bfs.webm -w 500 -h 400 --fps=5/1 --target-bitrate=2000 --good --threads=4".split())
   0


The resulting animation can be downloaded :download:`here
<bfs.webm>`, or played below if your browser supports WebM.

.. raw:: html

   <div style="text-align:center">
       <video id="sirs" src="../_downloads/bfs.webm" controls></video>
   </div>
