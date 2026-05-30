"use client";

import { useState, useEffect } from "react";
import { motion } from "framer-motion";
import { Leaf, Menu, X } from "lucide-react";
import { Button } from "@/components/ui/button";

const links = [
  { label: "Tjenester", href: "#tjenester" },
  { label: "Utforsk boligen", href: "#utforsk" },
  { label: "Om Forny", href: "#om-forny" },
];

export function Navbar() {
  const [scrolled, setScrolled] = useState(false);
  const [open, setOpen]         = useState(false);

  useEffect(() => {
    const handler = () => setScrolled(window.scrollY > 20);
    window.addEventListener("scroll", handler, { passive: true });
    return () => window.removeEventListener("scroll", handler);
  }, []);

  const scrollTo = (href: string) => {
    setOpen(false);
    document.querySelector(href)?.scrollIntoView({ behavior: "smooth" });
  };

  return (
    <motion.header
      initial={{ y: -60, opacity: 0 }}
      animate={{ y: 0, opacity: 1 }}
      transition={{ duration: 0.5 }}
      className={`fixed top-0 left-0 right-0 z-50 transition-all duration-300 ${
        scrolled ? "bg-white/95 backdrop-blur-md shadow-sm" : "bg-transparent"
      }`}
    >
      <div className="container mx-auto px-5 md:px-8 h-16 flex items-center justify-between">
        {/* Logo */}
        <button onClick={() => scrollTo("#topp")} className="flex items-center gap-2 group">
          <div className="w-8 h-8 rounded-lg bg-green-700 flex items-center justify-center">
            <Leaf size={16} className="text-white" />
          </div>
          <span className="text-xl font-bold text-slate-900 tracking-tight">forny</span>
        </button>

        {/* Desktop nav */}
        <nav className="hidden md:flex items-center gap-8">
          {links.map((l) => (
            <button
              key={l.href}
              onClick={() => scrollTo(l.href)}
              className="text-sm text-slate-600 hover:text-green-700 font-medium transition-colors"
            >
              {l.label}
            </button>
          ))}
        </nav>

        <div className="hidden md:block">
          <Button
            onClick={() => scrollTo("#kontakt")}
            className="bg-green-700 hover:bg-green-800 text-white text-sm font-semibold px-5"
          >
            Få gratis tilbud
          </Button>
        </div>

        {/* Mobile toggle */}
        <button className="md:hidden" onClick={() => setOpen(!open)} aria-label="Meny">
          {open ? <X size={22} /> : <Menu size={22} />}
        </button>
      </div>

      {/* Mobile menu */}
      {open && (
        <div className="md:hidden bg-white border-t border-slate-100 px-5 pb-5 pt-3 flex flex-col gap-4">
          {links.map((l) => (
            <button
              key={l.href}
              onClick={() => scrollTo(l.href)}
              className="text-left text-base text-slate-700 font-medium py-1"
            >
              {l.label}
            </button>
          ))}
          <Button
            onClick={() => scrollTo("#kontakt")}
            className="bg-green-700 hover:bg-green-800 text-white w-full mt-2"
          >
            Få gratis tilbud
          </Button>
        </div>
      )}
    </motion.header>
  );
}
