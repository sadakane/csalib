require 'mkmf'
$objs = %w{
  ruby.o ../csa.a
}
create_makefile('csa')
