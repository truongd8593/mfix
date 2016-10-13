# Common to fluid and solid phases

from __future__ import print_function, absolute_import, unicode_literals, division

class SpeciesHandler(object):
    def species_all_aliases(self):
        for (name, data) in self.fluid_species.items():
            data.get('alias', name)
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


    def find_species_phase(self, species_or_alias):
        if not self.fluid_solver_disabled:
            if species_or_alias in self.fluid_species:
                return 0
            for (name, data) in self.fluid_species.items():
                if data.get('alias') == species_or_alias:
                    return 0

        for (p, s) in self.solids_species.items():
            if species_or_alias in s:
                return p
            for (name, data) in s.items():
                if data.get('alias') == species_or_alias:
                    return p


    def species_mol_weight(self, species_or_alias):
        if not self.fluid_solver_disabled:
            for (name, data) in self.fluid_species.items():
                alias = data.get('alias', name)
                if alias == species_or_alias:
                    return data.get('mol_weight')
        for (p, subdict) in self.solids_species.items():
            for (name, data) in subdict.items():
                alias = data.get('alias', name)
                if alias == species_or_alias:
                    return data.get('mol_weight')

        return None


    def species_of_phase(self, p):
        if p == 0:
            return [data.get('alias', name)
                    for (name, data) in self.fluid_species.items()]
        else:
            return [data.get('alias', name)
                    for (name, data) in self.solids_species[p].items()]
