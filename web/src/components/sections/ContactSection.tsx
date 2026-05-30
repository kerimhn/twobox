"use client";

import { useState } from "react";
import { motion } from "framer-motion";
import { Send, CheckCircle2 } from "lucide-react";
import { Button } from "@/components/ui/button";

export function ContactSection() {
  const [sent, setSent] = useState(false);

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    setSent(true);
  };

  return (
    <section id="kontakt" className="py-24 bg-slate-950 relative overflow-hidden">
      {/* Background accent */}
      <div className="absolute inset-0 pointer-events-none">
        <div className="absolute top-0 left-1/2 -translate-x-1/2 w-[600px] h-[600px] bg-green-900/20 rounded-full blur-3xl" />
      </div>

      <div className="relative z-10 container mx-auto px-5 md:px-8">
        <div className="max-w-2xl mx-auto text-center mb-12">
          <span className="text-sm font-semibold uppercase tracking-widest text-green-400">Kom i gang</span>
          <h2 className="text-4xl font-bold text-white mt-2 mb-4">
            Få en gratis og uforpliktende prat
          </h2>
          <p className="text-slate-400 text-lg">
            Fortell oss litt om boligen din, så tar vi kontakt for å gå gjennom
            de beste energitiltakene for deg.
          </p>
        </div>

        <motion.div
          initial={{ opacity: 0, y: 24 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          className="max-w-lg mx-auto"
        >
          {sent ? (
            <div className="bg-green-900/40 border border-green-700/50 rounded-2xl p-10 text-center">
              <CheckCircle2 size={48} className="text-green-400 mx-auto mb-4" />
              <h3 className="text-xl font-bold text-white mb-2">Takk for henvendelsen!</h3>
              <p className="text-slate-400">
                Vi tar kontakt innen 1–2 virkedager for en uforpliktende prat.
              </p>
            </div>
          ) : (
            <form
              onSubmit={handleSubmit}
              className="bg-white/5 backdrop-blur-md border border-white/10 rounded-2xl p-8 space-y-4"
            >
              <div className="grid sm:grid-cols-2 gap-4">
                <Field label="Fornavn" name="firstname" required />
                <Field label="Etternavn" name="lastname" required />
              </div>
              <Field label="E-post" name="email" type="email" required />
              <Field label="Telefon (valgfritt)" name="phone" type="tel" />
              <div>
                <label className="block text-sm text-slate-300 font-medium mb-1.5">
                  Hva ønsker du å gjøre?
                </label>
                <textarea
                  name="message"
                  rows={4}
                  placeholder="Fortell oss om boligen og hva du tenker på…"
                  className="w-full rounded-xl bg-white/10 border border-white/15 text-white placeholder-slate-500 px-4 py-3 text-sm focus:outline-none focus:border-green-500 transition resize-none"
                />
              </div>
              <Button
                type="submit"
                className="w-full bg-green-600 hover:bg-green-500 text-white font-semibold text-base py-6 rounded-xl gap-2"
              >
                <Send size={16} /> Send forespørsel
              </Button>
              <p className="text-center text-xs text-slate-500">
                Helt uforpliktende. Vi deler ikke informasjonen din.
              </p>
            </form>
          )}
        </motion.div>
      </div>
    </section>
  );
}

function Field({
  label, name, type = "text", required = false,
}: {
  label: string; name: string; type?: string; required?: boolean;
}) {
  return (
    <div>
      <label htmlFor={name} className="block text-sm text-slate-300 font-medium mb-1.5">
        {label}
      </label>
      <input
        id={name}
        name={name}
        type={type}
        required={required}
        className="w-full rounded-xl bg-white/10 border border-white/15 text-white placeholder-slate-500 px-4 py-3 text-sm focus:outline-none focus:border-green-500 transition"
        placeholder={label}
      />
    </div>
  );
}
