"use client";

import { motion } from "framer-motion";
import { ArrowDown, Leaf } from "lucide-react";
import { Button } from "@/components/ui/button";

export function HeroSection() {
  const scroll = (id: string) =>
    document.getElementById(id)?.scrollIntoView({ behavior: "smooth" });

  return (
    <section
      id="topp"
      className="relative min-h-screen w-full flex items-center justify-center overflow-hidden bg-slate-950"
    >
      {/* Animated SVG mesh background */}
      <div className="absolute inset-0 pointer-events-none opacity-20">
        <svg className="w-full h-full" viewBox="0 0 1200 800" preserveAspectRatio="xMidYMid slice">
          {Array.from({ length: 20 }, (_, i) => (
            <motion.line
              key={i}
              x1={i * 60}
              y1={0}
              x2={i * 60 + 200}
              y2={800}
              stroke="#22c55e"
              strokeWidth="0.5"
              initial={{ opacity: 0.1 }}
              animate={{ opacity: [0.1, 0.4, 0.1] }}
              transition={{ duration: 4 + i * 0.3, repeat: Infinity, delay: i * 0.2 }}
            />
          ))}
          {Array.from({ length: 12 }, (_, i) => (
            <motion.circle
              key={`c${i}`}
              cx={100 + i * 100}
              cy={100 + (i % 4) * 180}
              r={40 + i * 8}
              stroke="#16a34a"
              strokeWidth="0.5"
              fill="none"
              initial={{ scale: 0.8, opacity: 0.2 }}
              animate={{ scale: [0.8, 1.2, 0.8], opacity: [0.2, 0.5, 0.2] }}
              transition={{ duration: 6 + i * 0.5, repeat: Infinity, delay: i * 0.4 }}
            />
          ))}
        </svg>
      </div>

      {/* Gradient overlay */}
      <div className="absolute inset-0 bg-gradient-to-br from-slate-950 via-slate-900 to-green-950 opacity-90" />

      {/* Content */}
      <div className="relative z-10 container mx-auto px-5 md:px-8 text-center">
        <motion.div
          initial={{ opacity: 0, y: 30 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.8, delay: 0.2 }}
          className="max-w-3xl mx-auto"
        >
          {/* Badge */}
          <motion.div
            initial={{ opacity: 0, scale: 0.9 }}
            animate={{ opacity: 1, scale: 1 }}
            transition={{ duration: 0.5 }}
            className="inline-flex items-center gap-2 bg-green-900/60 border border-green-700/50 rounded-full px-4 py-1.5 mb-8"
          >
            <Leaf size={13} className="text-green-400" />
            <span className="text-sm text-green-300 font-medium">Energieffektivisering for norske boliger</span>
          </motion.div>

          {/* Headline */}
          <h1 className="text-5xl sm:text-6xl md:text-7xl font-bold mb-6 text-white tracking-tight leading-none">
            {"Forny din bolig.".split(" ").map((word, wi) => (
              <span key={wi} className="inline-block mr-4 last:mr-0">
                {word.split("").map((letter, li) => (
                  <motion.span
                    key={`${wi}-${li}`}
                    initial={{ y: 60, opacity: 0 }}
                    animate={{ y: 0, opacity: 1 }}
                    transition={{
                      delay: 0.4 + wi * 0.1 + li * 0.025,
                      type: "spring",
                      stiffness: 160,
                      damping: 22,
                    }}
                    className="inline-block"
                  >
                    {letter}
                  </motion.span>
                ))}
              </span>
            ))}
            <br />
            <span className="text-green-400">Spar penger.</span>
          </h1>

          {/* Subtext */}
          <motion.p
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 1, duration: 0.7 }}
            className="text-lg md:text-xl text-slate-300 mb-10 max-w-xl mx-auto leading-relaxed"
          >
            Vi finner de mest lønnsomme energitiltakene for din bolig og kobler
            deg med riktig håndverker. Ingen oppgave er for liten.
          </motion.p>

          {/* CTAs */}
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 1.2, duration: 0.6 }}
            className="flex flex-col sm:flex-row gap-4 justify-center"
          >
            <Button
              onClick={() => scroll("utforsk")}
              className="bg-green-600 hover:bg-green-500 text-white text-base font-semibold px-8 py-6 rounded-xl"
            >
              Utforsk din bolig
            </Button>
            <Button
              onClick={() => scroll("kontakt")}
              variant="outline"
              className="border-slate-600 text-slate-200 hover:bg-slate-800 hover:text-white text-base font-semibold px-8 py-6 rounded-xl"
            >
              Kontakt oss gratis
            </Button>
          </motion.div>

          {/* Stats */}
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={{ delay: 1.6 }}
            className="mt-16 grid grid-cols-3 gap-6 max-w-lg mx-auto text-center"
          >
            {[
              { value: "50+", label: "Tiltak ranket" },
              { value: "30%", label: "Snittbesparelse" },
              { value: "100%", label: "Uforpliktende" },
            ].map((s) => (
              <div key={s.label}>
                <div className="text-2xl font-bold text-green-400">{s.value}</div>
                <div className="text-xs text-slate-400 mt-0.5">{s.label}</div>
              </div>
            ))}
          </motion.div>
        </motion.div>
      </div>

      {/* Scroll indicator */}
      <motion.div
        animate={{ y: [0, 8, 0] }}
        transition={{ repeat: Infinity, duration: 2 }}
        className="absolute bottom-8 left-1/2 -translate-x-1/2 text-slate-500"
      >
        <ArrowDown size={20} />
      </motion.div>
    </section>
  );
}
