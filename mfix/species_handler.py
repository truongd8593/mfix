# Common to fluid and solid phases

from __future__ import print_function, absolute_import, unicode_literals, division

class SpeciesHandler(object):
    def species_all_aliases(self):
        for k in self.fluid_species.keys():
            yield k
        for phase in self.solids_species.values():
            for k in phase.keys():
                yield k

    def species_alias_unique(self, alias):
        return alias not in self.species_all_aliases()

    def species_make_alias(self, species):
        if self.species_alias_unique(species):
            return species

        # Strip _nn suffix
        i = species.rindex('_')
        count = 0
        if i > 0 and species[i+1:].isdigit():
            count = 1 + int(species[i+1:])
            species = species[:i]

        while True:
            alias = '%s_%s' % (species, count) # Prefer underscore so we don't get things like H2O2
            if self.species_alias_unique(alias):
                return alias
            count += 1
