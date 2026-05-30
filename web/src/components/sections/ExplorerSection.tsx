"use client";

import { useCallback, useState } from "react";
import dynamic from "next/dynamic";
import { UpgradePanel } from "@/components/house/UpgradePanel";
import { UPGRADE_LIST } from "@/lib/energy-upgrades";
import { motion } from "framer-motion";

const HouseScene = dynamic(() => import("@/components/house/HouseScene"), {
  ssr: false,
  loading: () => (
    <div className="w-full h-full flex items-center justify-center bg-[#F0FBF4]">
      <div className="text-center">
        <div className="w-10 h-10 border-2 border-green-600 border-t-transparent rounded-full animate-spin mx-auto mb-3" />
        <p className="text-slate-400 text-sm">Laster 3D-modell…</p>
      </div>
    </div>
  ),
});

const PART_LABELS: Record<string, string> = {
  roof:       "Tak",
  walls:      "Vegger",
  windows:    "Vinduer",
  door:       "Ytterdør",
  foundation: "Gulv",
  solar:      "Solceller",
  heatpump:   "Varmepumpe",
};

export function ExplorerSection() {
  const [selected, setSelected] = useState<string | null>(null);
  const [hovered,  setHovered]  = useState<string | null>(null);

  const handleClose = useCallback(() => setSelected(null), []);

  return (
    <section id="utforsk" className="py-20 bg-slate-50">
      <div className="container mx-auto px-5 md:px-8 text-center mb-10">
        <span className="text-sm font-semibold uppercase tracking-widest text-green-700">
          Interaktiv bolig
        </span>
        <h2 className="text-4xl font-bold text-slate-900 mt-2 mb-4">
          Klikk på boligen – se hva du kan spare
        </h2>
        <p className="text-slate-500 max-w-xl mx-auto">
          Roter, zoom og klikk på de ulike delene av huset for å se hvilke
          energitiltak Forny tilbyr – med kostnader og besparelser.
        </p>
      </div>

      <div className="max-w-5xl mx-auto px-5 md:px-8">
        {/* Part pill shortcuts */}
        <div className="flex flex-wrap gap-2 justify-center mb-6">
          {Object.entries(PART_LABELS).map(([id, label]) => (
            <button
              key={id}
              onClick={() => setSelected(selected === id ? null : id)}
              className={`px-3 py-1.5 rounded-full text-xs font-semibold border transition-all duration-150 ${
                selected === id
                  ? "bg-green-700 text-white border-green-700"
                  : "bg-white text-slate-600 border-slate-200 hover:border-green-400 hover:text-green-700"
              }`}
            >
              {label}
            </button>
          ))}
        </div>

        {/* Canvas container */}
        <div className="relative h-[560px] rounded-3xl overflow-hidden shadow-2xl border border-slate-200">
          <HouseScene onPartSelect={setSelected} onPartHover={setHovered} />
          <UpgradePanel selected={selected} hovered={hovered} onClose={handleClose} />

          {/* Controls hint */}
          <div className="absolute bottom-4 left-4 flex items-center gap-3 bg-white/80 backdrop-blur-sm rounded-xl px-3 py-2 text-xs text-slate-500 select-none">
            <span>🖱 Dra for å rotere</span>
            <span>⚲ Scroll for å zoome</span>
            <span>👆 Klikk for å utforske</span>
          </div>
        </div>

        {/* Upgrade summary cards */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6 }}
          className="mt-8 grid grid-cols-2 md:grid-cols-4 gap-3"
        >
          {UPGRADE_LIST.slice(0, 4).map((u) => (
            <button
              key={u.id}
              onClick={() => setSelected(selected === u.id ? null : u.id)}
              className={`text-left p-4 rounded-2xl border transition-all duration-200 ${
                selected === u.id
                  ? "border-green-500 bg-green-50 shadow-md"
                  : "border-slate-200 bg-white hover:border-green-300 hover:shadow"
              }`}
            >
              <div
                className="w-3 h-3 rounded-full mb-2"
                style={{ backgroundColor: u.highlightColor }}
              />
              <p className="text-xs font-bold text-slate-700 leading-tight">{u.name}</p>
              <p className="text-xs text-green-700 font-semibold mt-1">{u.savings}</p>
            </button>
          ))}
        </motion.div>
      </div>
    </section>
  );
}
