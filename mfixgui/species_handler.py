# Common to fluid and solid phases

from __future__ import print_function, absolute_import, unicode_literals, division

class SpeciesHandler(object):
    def species_all_aliases(self):
        for (name, data) in self.fluid_species.items():
            yield data.get('alias', name)
        for phase in self.solids_species.values():
            for (name, data) in phase.items():
                yield data.get('alias', name)

    def species_alias_unique(self, alias):
        return alias not in self.species_all_aliases()

    def species_make_alias(self, species):
        if self.species_alias_unique(species):
            return species
        count = 0
        # Strip _nn suffix
        if '_' in species:
            i = species.rindex('_')
            if i > 0 and species[i+1:].isdigit():
                count = 1 + int(species[i+1:])
                species = species[:i]

        while True:
            alias = '%s_%s' % (species, count) # Prefer underscore so we don't get things like H2O2
            if self.species_alias_unique(alias):
                return alias
            count += 1


    def find_species_phase(self, species):
        if species in self.fluid_species:
            return 0

        for (p, s) in self.solids_species.items():
            if species in s:
                return p


    def species_mol_weight(self, species):
        if species in self.fluid_species:
            return self.fluid_species[species].get('mol_weight')
        for (p, subdict) in self.solids_species.items():
            if species in subdict:
                return subdict[species].get('mol_weight')

        return None


    def species_of_phase(self, p):
        if p is None:
            return []
        elif p == 0:
            return self.fluid_species.keys()
        else:
            return self.solids_species[p].keys()
