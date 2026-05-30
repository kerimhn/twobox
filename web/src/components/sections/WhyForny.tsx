"use client";

import { motion } from "framer-motion";
import { Award, BarChart2, Leaf, PhoneCall } from "lucide-react";

const reasons = [
  {
    icon: BarChart2,
    title: "Rangert etter lønnsomhet",
    description:
      "Vi analyserer kostnaden og besparelsen for hvert tiltak, og presenterer dem i rekkefølge – mest lønnsomt øverst.",
    color: "text-blue-600",
    bg: "bg-blue-50",
  },
  {
    icon: Award,
    title: "Sertifiserte byggmestre",
    description:
      "Vi samarbeider kun med godkjente håndverkere med dokumentert erfaring innen energieffektivisering.",
    color: "text-amber-600",
    bg: "bg-amber-50",
  },
  {
    icon: Leaf,
    title: "Enova-støtte inkludert",
    description:
      "Vi hjelper deg med å søke Enova-støtte slik at du får mest mulig igjen for investeringen din.",
    color: "text-green-700",
    bg: "bg-green-50",
  },
  {
    icon: PhoneCall,
    title: "Ingen oppgave for liten",
    description:
      "Enten det er ett vindu eller en total energioppgradering – vi tar deg seriøst og finner den beste løsningen.",
    color: "text-purple-600",
    bg: "bg-purple-50",
  },
];

export function WhyForny() {
  return (
    <section id="om-forny" className="py-24 bg-white">
      <div className="container mx-auto px-5 md:px-8">
        <div className="text-center mb-16">
          <span className="text-sm font-semibold uppercase tracking-widest text-green-700">Hvorfor Forny</span>
          <h2 className="text-4xl font-bold text-slate-900 mt-2 mb-4">
            Lønnsomhet, ikke bare miljøtanker
          </h2>
          <p className="text-slate-500 max-w-xl mx-auto text-lg">
            Vi hjelper deg å ta smarte beslutninger som sparer deg penger – og tilfeldigvis er bra for klimaet.
          </p>
        </div>

        <div className="grid sm:grid-cols-2 lg:grid-cols-4 gap-6 max-w-5xl mx-auto">
          {reasons.map((r, i) => (
            <motion.div
              key={i}
              initial={{ opacity: 0, y: 24 }}
              whileInView={{ opacity: 1, y: 0 }}
              viewport={{ once: true }}
              transition={{ delay: i * 0.1, duration: 0.5 }}
              className="p-6 rounded-2xl border border-slate-100 hover:border-green-200 hover:shadow-md transition-all duration-200"
            >
              <div className={`w-12 h-12 ${r.bg} rounded-xl flex items-center justify-center mb-4`}>
                <r.icon size={22} className={r.color} />
              </div>
              <h3 className="text-base font-bold text-slate-900 mb-2">{r.title}</h3>
              <p className="text-sm text-slate-500 leading-relaxed">{r.description}</p>
            </motion.div>
          ))}
        </div>

        {/* Trust strip */}
        <div className="mt-16 py-8 border-t border-slate-100 flex flex-col md:flex-row items-center justify-center gap-10 text-center text-sm text-slate-400">
          {[
            "Gratis og uforpliktende",
            "Godkjent av Enova",
            "Lokal håndverker",
            "Hele landet",
          ].map((item) => (
            <div key={item} className="flex items-center gap-2">
              <div className="w-1.5 h-1.5 rounded-full bg-green-500" />
              <span>{item}</span>
            </div>
          ))}
        </div>
      </div>
    </section>
  );
}
