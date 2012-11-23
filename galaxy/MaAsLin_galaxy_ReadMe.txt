Installation instructions for maaslin in a galaxy environment.
These instructions require the Mercurial versioning system, galaxy, and an internet connection.

1. In the  "galaxy-dist/tools" directory install maaslin by typing in a terminal:
hg clone https://bitbucket.org/chuttenh/maaslin

2. Update member tool_conf.xml  in the galaxy directory adding the following: 
  <section name="maaslin" id="maaslin">
    <tool file="maaslin/galaxy/maaslin_input.xml"/>
    <tool file="maaslin/galaxy/maaslin.xml"/>
  </section>

3. Update member datatypes_conf.xml  in the galaxy directory adding the following:
	<datatype extension="maaslin" type="galaxy.datatypes.data:Text" subclass="true" display_in_upload="true"/>

4. Copy members Figure1-Overview.png and Maaslin_Output.png  to /galaxy/static/images 
5. Recycle galaxy

