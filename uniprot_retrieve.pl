#!/usr/bin/perl


#program from http://www.uniprot.org/help/programmatic_access , MU Sun Apr 19 13:44:32 BST 2015
use strict;
use warnings;
use LWP::UserAgent;

my $list = $ARGV[0]; # File containg list of UniProt identifiers.

my $base = 'http://www.uniprot.org';
my $tool = 'batch';

my $contact = 'Martin.Urban@rothamsted.ac.uk'; # Please set your email address here to help us debug in case of problems.
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
push @{$agent->requests_redirectable}, 'POST';

my $response = $agent->post("$base/$tool/",
                            [ 'file' => [$list],
                              'format' => 'fasta',
                            ],
                            'Content_Type' => 'form-data');

while (my $wait = $response->header('Retry-After')) {
  print STDERR "Waiting ($wait)...\n";
  sleep $wait;
  $response = $agent->get($response->base);
}

$response->is_success ?
  print $response->content :
  die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";
