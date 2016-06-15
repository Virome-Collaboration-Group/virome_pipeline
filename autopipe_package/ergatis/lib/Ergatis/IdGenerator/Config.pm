package Ergatis::IdGenerator::Config;

## This config file defines the module that will be used for id generation by ergatis
## The specified class should at least implement all methods defined in Ergatis::IdGenerator

use Ergatis::IdGenerator::RandomIdGenerator;
our $class = 'Ergatis::IdGenerator::RandomIdGenerator';

1;
