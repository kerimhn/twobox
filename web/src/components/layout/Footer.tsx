import { Leaf } from "lucide-react";

export function Footer() {
  return (
    <footer className="bg-slate-900 text-slate-400 py-12">
      <div className="container mx-auto px-5 md:px-8">
        <div className="grid md:grid-cols-3 gap-8 mb-10">
          {/* Brand */}
          <div>
            <div className="flex items-center gap-2 mb-3">
              <div className="w-7 h-7 rounded-md bg-green-600 flex items-center justify-center">
                <Leaf size={14} className="text-white" />
              </div>
              <span className="text-lg font-bold text-white">forny</span>
            </div>
            <p className="text-sm leading-relaxed">
              Vi hjelper norske huseiere med å finne de mest lønnsomme
              energitiltakene – ingen oppgave er for liten.
            </p>
          </div>

          {/* Links */}
          <div>
            <h4 className="text-sm font-semibold text-white mb-3 uppercase tracking-wider">Tjenester</h4>
            <ul className="space-y-2 text-sm">
              {["Loftsisolasjon", "Veggsisolasjon", "Vindusskifte", "Varmepumpe", "Solcelleanlegg"].map((s) => (
                <li key={s}>{s}</li>
              ))}
            </ul>
          </div>

          {/* Contact */}
          <div>
            <h4 className="text-sm font-semibold text-white mb-3 uppercase tracking-wider">Kontakt</h4>
            <ul className="space-y-2 text-sm">
              <li>hei@forny.org</li>
              <li>forny.org</li>
              <li className="pt-1 text-xs">
                Vi samarbeider med sertifiserte byggmestre over hele landet.
              </li>
            </ul>
          </div>
        </div>

        <div className="border-t border-slate-800 pt-6 flex flex-col md:flex-row justify-between items-center gap-3 text-xs">
          <p>© {new Date().getFullYear()} Forny AS. Alle rettigheter forbeholdt.</p>
          <p>Designet for et mer energieffektivt Norge.</p>
        </div>
      </div>
    </footer>
  );
}
