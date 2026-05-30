"use client";

import { AnimatePresence, motion } from "framer-motion";
import { X, CheckCircle2, Zap, Clock, Leaf } from "lucide-react";
import { Button } from "@/components/ui/button";
import { UPGRADES } from "@/lib/energy-upgrades";

interface UpgradePanelProps {
  selected: string | null;
  hovered:  string | null;
  onClose:  () => void;
}

export function UpgradePanel({ selected, hovered, onClose }: UpgradePanelProps) {
  const activeId = selected ?? hovered;
  const upgrade  = activeId ? UPGRADES[activeId] : null;
  const isPinned = !!selected;

  return (
    <AnimatePresence mode="wait">
      {upgrade && (
        <motion.div
          key={activeId}
          initial={{ opacity: 0, x: 32, scale: 0.97 }}
          animate={{ opacity: 1, x: 0,  scale: 1 }}
          exit={{    opacity: 0, x: 32, scale: 0.97 }}
          transition={{ type: "spring", stiffness: 320, damping: 28 }}
          className="absolute top-4 right-4 w-72 rounded-2xl overflow-hidden shadow-2xl border border-white/60 bg-white/95 backdrop-blur-md"
        >
          {/* Colour accent bar */}
          <div className="h-1" style={{ backgroundColor: upgrade.highlightColor }} />

          <div className="p-5">
            {isPinned && (
              <button
                onClick={onClose}
                className="absolute top-3 right-3 text-slate-400 hover:text-slate-600 transition-colors"
                aria-label="Lukk"
              >
                <X size={17} />
              </button>
            )}

            <p className="text-[11px] font-semibold tracking-widest uppercase text-slate-400 mb-0.5">
              {upgrade.subtitle}
            </p>
            <h3 className="text-lg font-bold text-slate-900 leading-tight">
              {upgrade.name}
            </h3>

            <p className="mt-2 text-sm text-slate-600 leading-relaxed">
              {upgrade.description}
            </p>

            {isPinned && (
              <>
                <div className="mt-4 space-y-2">
                  <StatRow icon={<Zap size={13} className="text-amber-500" />} label="Besparelse" value={upgrade.savings} />
                  <StatRow icon={<CheckCircle2 size={13} className="text-green-600" />} label="Kostnad" value={upgrade.costRange} />
                  <StatRow icon={<Clock size={13} className="text-blue-500" />} label="Tilbakebetaling" value={upgrade.payback} />
                  <StatRow icon={<Leaf size={13} className="text-emerald-500" />} label="Enova-støtte" value={upgrade.enovaSupport} highlight />
                </div>

                <div className="mt-5">
                  <Button
                    className="w-full bg-green-700 hover:bg-green-800 text-white text-sm font-semibold"
                    onClick={() => {
                      document.getElementById("kontakt")?.scrollIntoView({ behavior: "smooth" });
                      onClose();
                    }}
                  >
                    Få gratis tilbud
                  </Button>
                </div>
              </>
            )}

            {!isPinned && (
              <p className="mt-3 text-xs text-slate-400">Klikk for å se detaljer og pris</p>
            )}
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}

function StatRow({
  icon, label, value, highlight = false,
}: {
  icon: React.ReactNode; label: string; value: string; highlight?: boolean;
}) {
  return (
    <div className={`flex items-center gap-2 text-sm rounded-lg px-3 py-2 ${highlight ? "bg-emerald-50" : "bg-slate-50"}`}>
      {icon}
      <span className="text-slate-500 flex-1">{label}</span>
      <span className="font-semibold text-slate-800">{value}</span>
    </div>
  );
}
