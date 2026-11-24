// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/dst/dst_writer.hpp"

#include <iomanip>

namespace delphipp::dst {

void write_csv_header(std::ostream& out) {
  out << "run,event,particle_id,charge,pid,px,py,pz,energy,"
         "impact_rphi_cm,impact_z_cm,ecal_gev,hcal_gev,tpc_dedx,muon_tag\n";
}

void write_csv_event(std::ostream& out, const DstEvent& ev) {
  out << std::setprecision(6);
  for (const auto& p : ev.particles) {
    out << ev.run << ',' << ev.event << ',' << p.id << ',' << p.charge << ','
        << static_cast<std::int32_t>(p.pid) << ',' << p.px << ',' << p.py << ','
        << p.pz << ',' << p.energy << ',' << p.impact_rphi_cm << ','
        << p.impact_z_cm << ',' << p.ecal_energy_gev << ',' << p.hcal_energy_gev
        << ',' << p.tpc_dedx << ',' << (p.muon_tag ? 1 : 0) << '\n';
  }
}

void write_ndjson_event(std::ostream& out, const DstEvent& ev) {
  out << ev.to_json() << '\n';
}

void write_hepmc3_header(std::ostream& out) {
  out << "HepMC::Version 3.02.06\n";
  out << "HepMC::Asciiv3-START_EVENT_LISTING\n";
}

void write_hepmc3_event(std::ostream& out, const DstEvent& ev) {
  out << std::setprecision(6);
  // E event_number num_vertices num_particles
  out << "E " << ev.event << " 1 " << ev.particles.size() << '\n';
  out << "U GEV CM\n";
  // V -1 0 [x, y, z, t]
  out << "V -1 0 [" << ev.primary_vertex.xyz_cm[0] << ','
      << ev.primary_vertex.xyz_cm[1] << ','
      << ev.primary_vertex.xyz_cm[2] << ",0]\n";
  for (const auto& p : ev.particles) {
    // P id parent_vtx pid px py pz E m status
    const std::int32_t pdg =
        (p.charge < 0 && (p.pid == ParticleTag::Electron || p.pid == ParticleTag::Muon))
            ? static_cast<std::int32_t>(p.pid)
            : (p.charge > 0 && (p.pid == ParticleTag::Electron || p.pid == ParticleTag::Muon))
                  ? -static_cast<std::int32_t>(p.pid)
                  : ((p.charge < 0) ? -static_cast<std::int32_t>(p.pid)
                                    : static_cast<std::int32_t>(p.pid));
    out << "P " << p.id << " -1 " << pdg << ' ' << p.px << ' ' << p.py << ' '
        << p.pz << ' ' << p.energy << " 0 1\n";
  }
}

void write_hepmc3_footer(std::ostream& out) {
  out << "HepMC::Asciiv3-END_EVENT_LISTING\n";
}

}  // namespace delphipp::dst
