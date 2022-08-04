var TidyTree = (function () {
  'use strict';

  (function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
    typeof define === 'function' && define.amd ? define(['exports'], factory) :
    (global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.patristic = {}));
  }(undefined, (function (exports) {
    /**
     * The [SemVer](https://semver.org/) version string of the patristic library
     * @type {String} A string specifying the current version of the Patristic Library.
     * If not given, the version of patristic you are using if less than or equal to 0.2.2.
     * @example
     * console.log(patristic.version);
     */
    const version = "0.5.7";

    /**
     * A class for representing Branches in trees.
     * It's written predominantly for phylogenetic trees (hence the
     * [Newick parser](#parseNewick),
     * [neighbor-joining implementation](#parseMatrix), etc.), but could
     * conceivably be useful for representing other types of trees as well.
     * @param {Object} [data] An object containing data you wish to assign to
     * this Branch object. In particular, intended to overwrite the default
     * attributes of a Branch, namely `id`, `parent`, `length`, and `children`.
     * @constructor
     */
    function Branch(data, children) {
      if (!data) data = {};
      if (!children) children = d => d.children;
      Object.assign(this, {
        _guid: guid(),
        id: data.id || "",
        data: data,
        depth: data.depth || 0,
        height: data.height || 0,
        length: data.length || 0,
        parent: data.parent || null,
        children: children(data) || [],
        value: data.value || 1,
        respresenting: 1
      });
    }

    function guid(){
      return ([1e7] + -1e3 + -4e3 + -8e3 + -1e11).replace(/[018]/g, a => {
        return (a ^ ((Math.random() * 16) >> (a / 4))).toString(16);
      });
    }

    /**
     * Adds a new child to this Branch
     * @param  {(Branch|Object)} [data={}] The new Branch, or data to attach to it.
     * @return {Branch} The (possibly new) child Branch
     */
    Branch.prototype.addChild = function(data) {
      let c;
      if (data instanceof Branch) {
        c = data;
        c.parent = this;
      } else {
        if (!data) data = {};
        c = new Branch(
          Object.assign(data, {
            parent: this
          })
        );
      }
      this.children.push(c);
      return c;
    };

    /**
     * Adds a new parent to this Branch. This is a bit esoteric and generally not
     * recommended.
     * @param  {(Branch|Object)} [data={}] A Branch object, or the data to attach to one
     * @param  {Array} [siblings=[]] An array of Branches to be the children of the new parent Branch (i.e. siblings of this Branch)
     * @return {Branch} The Branch on which this was called
     */
    Branch.prototype.addParent = function(data, siblings) {
      if (!siblings) siblings = [];
      let c;
      if (data instanceof Branch) {
        c = data;
      } else {
        if (!data) data = {};
        c = new Branch(Object.assign(data));
      }
      siblings.forEach(sib => sib.setParent(c));
      c.children = [this].concat(siblings);
      this.parent = c;
      return this;
    };

    /**
     * Returns an array of Branches from this Branch to the root.
     * [d3-hierarchy compatibility method.](https://github.com/d3/d3-hierarchy#node_ancestors)
     * @type {Array} An array of Branches
     */
    Branch.prototype.ancestors = function() {
      return this.getAncestors(true);
    };

    /**
     * Returns a deep clone of the Branch on which it is called. Note that this does
     * not clone all descendants, rather than providing references to the existing
     * descendant Branches.
     * @return {Branch} A clone of the Branch on which it is called.
     */
    Branch.prototype.clone = function() {
      return parseJSON(this.toObject());
    };

    /**
     * All descendant Branches with near-zero length are excised
     * @return {Branch} The Branch on which this method was called.
     */
    Branch.prototype.consolidate = function() {
      return this.eachAfter(branch => {
        if (branch.isRoot() || branch.length >= 0.0005) return;
        if(branch.parent.id == ""){
          branch.parent.id = branch.id;
        } else {
          branch.parent.id += '+'+branch.id;
        }
        branch.excise();
      }).fixDistances();
    };

    /**
     * Returns a clone of the Branch on which it is called. Note that this also
     * clones all descendants, rather than providing references to the existing
     * descendant Branches. (For a shallow clone, see [Branch.clone](#clone).
     * Finally, the cloned Branch will become the root of the cloned tree, having a
     * parent of `null`.
     * [d3-hierarchy compatibility method.](https://github.com/d3/d3-hierarchy#node_copy)
     * @return {Branch} A clone of the Branch on which it is called.
     */
    Branch.prototype.copy = function() {
      let newThis = parseJSON(JSON.stringify(this));
      newThis.parent = null;
      return newThis.fixDistances();
    };

    /**
     * Sets the values of all nodes to be equal to the number of their descendants.
     * @return {Branch} The Branch on which it was called
     */
    Branch.prototype.count = function() {
      return this.sum(() => 1);
    };

    /**
     * Returns an array pf descendants, starting with this Branch.
     * [d3-hierarchy compatibility method.](https://github.com/d3/d3-hierarchy#node_descendants)
     * @type {Array} An Array of Branches, starting with this one.
     */
    Branch.prototype.descendants = function() {
      return this.getDescendants(true);
    };

    /**
     * Returns the depth of a given child, relative to the Branch on which it is
     * called.
     * @param  {(Branch|String)} descendant A descendant Branch (or `id` string
     * thereof)
     * @return {Number} The sum of the lengths of all Branches between the Branch on
     * which it is called and `descendant`. Throws an error if `descendant` is not a
     * descendant of this Branch.
     */
    Branch.prototype.depthOf = function(descendant) {
      let distance = 0;
      if (typeof descendant == "string")
        descendant = this.getDescendant(descendant);
      if (typeof descendant == "undefined")
        throw Error("Cannot compute depth of undefined descendant!");
      let current = descendant;
      while (current != this) {
        distance += current.length;
        current = current.parent;
      }
      return distance;
    };

    /**
     * Computes the patristic distance between `descendantA` and `descendantB`.
     * @param  {Branch} descendantA The Branch from which you wish to compute
     * distance
     * @param  {Branch} descendantB The Branch to which you wish to compute distance
     * @return {number} The patristic distance between the given descendants.
     */
    Branch.prototype.distanceBetween = function(descendantA, descendantB) {
      let mrca = descendantA.getMRCA(descendantB);
      return mrca.depthOf(descendantA) + mrca.depthOf(descendantB);
    };

    /**
     * Computes the patristic distance between `cousin` and the Branch on which
     * this method is called.
     * @param  {Branch} cousin The Branch to which you wish to compute distance
     * @return {number} The patristic distance between `cousin` and the Branch on
     * this method is called.
     */
    Branch.prototype.distanceTo = function(cousin) {
      let mrca = this.getMRCA(cousin);
      return mrca.depthOf(this) + mrca.depthOf(cousin);
    };

    /**
     * Visits each Branch descended from the Branch on which it is called in
     * [Breadth First Search](https://en.wikipedia.org/wiki/Breadth-first_search)
     * order and returns the Branch on which it was called.
     * @param  {Function} callback The function to be run on each Branch
     * @return {Branch} The Branch on which it was called.
     */
    Branch.prototype.each = function(callback) {
      let branch = this,
        next = [branch],
        current;
      while (next.length) {
        current = next.reverse();
        next = [];
        while ((branch = current.pop())) {
          callback(branch);
          branch.eachChild(child => next.push(child));
        }
      }
      return this;
    };

    /**
     * Visits each Branch descended from the Branch on which it is called in
     * [post-traversal order](https://en.wikipedia.org/wiki/Tree_traversal#Post-order)
     * and returns the Branch on which it was called.
     * @param  {Function} callback Function to run on each Branch
     * @return {Branch} The Branch on which it was called
     */
    Branch.prototype.eachAfter = function(callback) {
      this.eachChild(child => child.eachAfter(callback));
      callback(this);
      return this;
    };

    /**
     * Visits each Branch descended from the Branch on which it is called in
     * [pre-traversal order](https://en.wikipedia.org/wiki/Tree_traversal#Pre-order)
     * and returns the Branch on which it was called.
     * @param  {Function} callback [description]
     * @return {[type]}            [description]
     */
    Branch.prototype.eachBefore = function(callback) {
      callback(this);
      this.eachChild(child => child.eachBefore(callback));
      return this;
    };

    /**
     * Runs a function on each child of the Branch on which it is called.
     * @param  {Function} callback The function to run on each child.
     * @return {Branch} The Branch on which it was called.
     */
    Branch.prototype.eachChild = function(callback) {
      this.children.forEach(callback);
      return this;
    };

    /**
     * Excises the Branch on which it is called and updates its parent and children.
     * @return {Branch} The parent of the excised Branch.
     */
    Branch.prototype.excise = function() {
      if (this.isRoot() && this.children.length > 1) {
        throw new Error("Cannot excise a root Branch with multiple children.");
      }
      this.eachChild(child => {
        child.length += this.length;
        child.parent = this.parent;
        if (!this.isRoot()) this.parent.children.push(child);
      });
      this.parent.children.splice(this.parent.children.indexOf(this), 1);
      this.parent.representing++;
      return this.parent;
    };

    /**
     * Sets the distance values (height and depth) for each Branch
     * @return {Branch} The Branch on which it is called.
     */
    Branch.prototype.fixDistances = function() {
      let maxdepth = 0,
        root = this.getRoot();
      root.depth = 0;
      this.eachBefore(d => {
        if (d.isRoot()) return;
        d.depth = d.parent.depth + 1;
        if (d.depth > maxdepth) maxdepth = d.depth;
      }).eachAfter(d => {
        d.height = maxdepth - d.depth;
        d.value = d.value + d.children.reduce((a, c) => a + c.value, 0);
      });
      return this;
    };

    /**
     * Repairs incorrect links by recurively confirming that children reference
     * their parents, and correcting those references if they do not.
     *
     * If you need to call this, something has messed up the state of your tree
     * and you should be concerned about that. Just FYI. ¯\\_(ツ)_/¯
     * @param  {Boolean} nonrecursive Should this just fix the children of the
     * Branch on which it is called, or all descendants?
     * @return {Branch} The Branch on which it was called.
     */
    Branch.prototype.fixParenthood = function(nonrecursive) {
      this.children.forEach(child => {
        if (!child.parent) child.parent = this;
        if (child.parent !== this) child.parent = this;
        if (!nonrecursive && child.children.length > 0) {
          child.fixParenthood();
        }
      });
      return this;
    };

    /**
     * Reverses the order of (all of) the descendants of the Branch.
     * @return {Branch} The Branch on which this was called.
     */
    Branch.prototype.flip = function() {
      return this.each(c => c.rotate());
    };

    /**
     * Returns an Array of all the ancestors of the Branch on which it is called.
     * Note that this does not include itself. For all ancestors and itself, see
     * [Branch.ancestors](#ancestors)
     * @param {Boolean} includeSelf Should the Branch on which this is called be
     * included in the results?
     * @return {Array} Every Ancestor of the Branch on which it was called.
     */
    Branch.prototype.getAncestors = function(includeSelf) {
      let ancestors = includeSelf ? [this] : [];
      let current = this;
      while ((current = current.parent)) ancestors.push(current);
      return ancestors;
    };

    /**
     * Given an `childID`, returns the child with that id (or `undefined` if no such
     * child is present).
     * @param  {String} childID the ID of the child to return.
     * @return {(Branch|undefined)} The desired child Branch, or `undefined` if the
     * child doesn't exist.
     */
    Branch.prototype.getChild = function(childID) {
      if (!typeof childID == "string") throw Error("childID is not a String!");
      return this.children.find(c => c.id === childID);
    };

    /**
     * Given an id string, returns the descendant Branch with that ID, or
     * `undefined` if it doesn't exist.
     * @param  {String} id The id string of the Branch to find
     * @return {(Branch|undefined)} The descendant Branch, or `undefined` if it
     * doesn't exist
     */
    Branch.prototype.getDescendant = function(id) {
      if (this.id === id) return this;
      let children = this.children,
        n = children.length;
      if (children) {
        for (let i = 0; i < n; i++) {
          let descendant = children[i].getDescendant(id);
          if (descendant) return descendant;
        }
      }
    };

    /**
     * Returns an array of all Branches which are descendants of this Branch
     * @param {Boolean} [includeSelf] Is this not the Branch on which the user
     * called the function? This is used internally and should be ignored.
     * @return {Array} An array of all Branches descended from this Branch
     */
    Branch.prototype.getDescendants = function(includeSelf) {
      let descendants = includeSelf ? [this] : [];
      if (!this.isLeaf()) {
        this.children.forEach(child => {
          child.getDescendants(true).forEach(d => descendants.push(d));
        });
      }
      return descendants;
    };

    /**
     * Returns an array of all leaves which are descendants of this Branch.
     * Alias of [getLeaves](#getLeaves) for people whose strong suit isn't spelling.
     * @return {Array} An array of all leaves descended from this Branch
     */
    Branch.prototype.getLeafs = function() {
      return this.getLeaves();
    };

    /**
     * Returns an array of all leaves which are descendants of this Branch
     * See also: [getLeafs](#getLeafs)
     * @return {Array} An array of all leaves descended from this Branch
     */
    Branch.prototype.getLeaves = function() {
      if (this.isLeaf()) {
        return [this];
      } else {
        let descendants = [];
        this.children.forEach(child => {
          child.getLeaves().forEach(d => descendants.push(d));
        });
        return descendants;
      }
    };

    /**
     * Traverses the tree upward until it finds the Most Recent Common Ancestor
     * (i.e. the first Branch for which both the Branch on which it was called and
     * `cousin` are descendants).
     * @return {Branch} The Most Recent Common Ancestor of both the Branch on
     * which it was called and the `cousin`.
     */
    Branch.prototype.getMRCA = function(cousin) {
      let mrca = this;
      while (!mrca.hasDescendant(cousin)) {
        if (mrca.isRoot())
          throw Error(
            "Branch and cousin do not appear to share a common ancestor!"
          );
        mrca = mrca.parent;
      }
      return mrca;
    };

    /**
     * Traverses the tree upward until it finds the root Branch, and returns the
     * root.
     * @return {Branch} The root Branch of the tree
     */
    Branch.prototype.getRoot = function() {
      let branch = this;
      while (!branch.isRoot()) branch = branch.parent;
      return branch;
    };

    /**
     * Determines if a given Branch (or ID) is a child of this Branch
     * @param  {(Branch|String)} child The Branch (or the id thereof) to check for
     * @return {Boolean}
     */
    Branch.prototype.hasChild = function(child) {
      if (child instanceof Branch) {
        return this.children.includes(child);
      } else if (typeof child === "string") {
        return this.children.some(c => c.id === child);
      }
      throw Error(
        `Unknown type of child (${typeof child}) passed to Branch.hasChild!`
      );
    };

    /**
     * Checks to see if `descendant` is a descendant of the Branch on which this
     * method is called.
     * @param  {(Branch|String)} descendant Either the descendant Branch or its'
     * `id`.
     * @return {Boolean} True if `descendant` is descended from the Branch from
     * which this is called, otherwise false.
     */
    Branch.prototype.hasDescendant = function(descendant) {
      let descendants = this.getDescendants();
      if (descendant instanceof Branch) {
        return descendants.some(d => d === descendant);
      } else if (typeof descendant === "string") {
        return descendants.some(d => d.id === descendant);
      }
      throw Error("Unknown type of descendant passed to Branch.hasDescendant!");
    };

    /**
     * Checks to see if a Branch has a descendant leaf.
     * @return {Boolean} True if leaf is both a leaf and a descendant of the
     * Branch on which this method is called, False otherwise.
     */
    Branch.prototype.hasLeaf = function(leaf) {
      let leaves = this.getleaves();
      if (leaf instanceof Branch) {
        return leaves.includes(leaf);
      } else if (typeof leaf === "string") {
        return leaves.some(d => d.id === leaf);
      }
      throw Error("Unknown type of leaf passed to Branch.hasLeaf.");
    };

    /**
     * Swaps the branch on which it is called with its parent. This method is
     * probably only useful as an internal component of [Branch.reroot](#reroot).
     * @return {Branch} The Branch object on which it was called.
     */
    Branch.prototype.invert = function() {
      let oldParent = this.parent;
      if (oldParent) {
        let temp = this.parent.length;
        this.parent.length = this.length;
        this.length = temp;
        this.parent = oldParent.parent;
        this.children.push(oldParent);
        oldParent.parent = this;
        oldParent.children.splice(oldParent.children.indexOf(this), 1);
      } else {
        throw Error("Cannot invert root node!");
      }
      return this;
    };

    /**
     * Returns whether the Branch on which it is called is a child of a given parent
     * (or parent ID).
     * @param  {(Branch|String)} parent A Branch (or ID thereof) to test for
     * paternity of this Branch.
     * @return {Boolean} True is `parent` is the parent of this Branch, false
     * otherwise.
     */
    Branch.prototype.isChildOf = function(parent) {
      if (parent instanceof Branch) return this.parent === parent;
      if (typeof parent === "string") return this.parent.id === parent;
      throw Error("Unknown parent type passed to Branch.isChildOf");
    };

    /**
     * Tests whether this and each descendant Branch holds correct links to both
     * its parent and its children.
     * @return {Boolean} True if consistent, otherwise false
     */
    Branch.prototype.isConsistent = function() {
      if (!this.isRoot()) {
        if (!this.parent.children.includes(this)) return false;
      }
      if (!this.isLeaf()) {
        if (this.children.some(c => c.parent !== this)) return false;
        return this.children.every(c => c.isConsistent());
      }
      return true;
    };

    /**
     * Returns whether a given Branch is an ancestor of the Branch on which this
     * method is called. Uses recursive tree-climbing.
     * @param  {Branch} ancestor The Branch to check for ancestorhood
     * @return {Boolean} If this Branch is descended from `ancestor`
     */
    Branch.prototype.isDescendantOf = function(ancestor) {
      if (!ancestor || !this.parent) return false;
      if (this.parent === ancestor || this.parent.id === ancestor) return true;
      return this.parent.isDescendantOf(ancestor);
    };

    /**
     * Returns a boolean indicating if this Branch is a leaf (i.e. has no
     * children).
     * @return {Boolean} True is this Branch is a leaf, otherwise false.
     */
    Branch.prototype.isLeaf = function() {
      return this.children.length === 0;
    };

    /**
     * Returns a boolean indicating whether or not this Branch is olate.
     *
     * ...Just kidding!
     *
     * Isolates a Branch and its subtree (i.e. removes everything above it, making
     * it the root Branch). Similar to [Branch.remove](#remove), only it returns
     * the Branch on which it is called.
     * @return {Branch} The Branch object on which it was called.
     */
    Branch.prototype.isolate = function() {
      let index = this.parent.children.indexOf(this);
      this.parent.children.splice(index, 1);
      this.setParent(null);
      return this;
    };

    /**
     * Returns a boolean indicating if this Branch is the root of a tree (i.e. has
     * no parents).
     * @return {Boolean} True if this Branch is the root, otherwise false.
     */
    Branch.prototype.isRoot = function() {
      return this.parent === null;
    };

    /**
     * Returns the array of leaf nodes in traversal order; leaves are nodes with no
     * children. Alias of [Branch.getLeaves](#getLeaves) `cuz spelling is hard.
     * @type {Array} An Array of Branches which are descended from this Branch and
     * have no children.
     */
    Branch.prototype.leafs = function() {
      return this.getLeaves();
    };

    /**
     * Returns the array of leaf nodes in traversal order; leaves are nodes with no
     * children. Alias of [Branch.getLeaves](#getLeaves).
     * [d3-hierarchy compatibility method.](https://github.com/d3/d3-hierarchy#node_leaves)
     * @type {Array} An Array of Branches which are descended from this Branch and
     * have no children.
     */
    Branch.prototype.leaves = function() {
      return this.getLeaves();
    };

    /**
     * Returns an Array of links, which are plain javascript objects containing a
     * `source` attribute (which is a reference to the parent Branch) and a `target`
     * attribute (which is a reference to the child Branch).
     * [d3-hierarchy compatibility method](https://github.com/d3/d3-hierarchy#node_links)
     * @return {Array} An array of plain Javascript objects
     */
    Branch.prototype.links = function() {
      let links = [];
      this.each(d => {
        if (d.isRoot()) return;
        links.push({
          source: d.parent,
          target: d
        });
      });
      return links;
    };

    /**
     * Normalizes this and all descendant Branches `value` attributes to between
     * `newmin` and `newmax`. Note that normalize can function as its own inverse
     * when passed an original range. For example:
     * @example tree.normalize().normalize(1, tree.getDescendants().length + 1);
     * @param  {Number} newmin The desired minimum value.
     * @param  {Number} newmax The desired maximum value.
     * @return {Branch} The Branch on which it was called.
     */
    Branch.prototype.normalize = function(newmin, newmax) {
      if (typeof newmax !== "number") newmax = 1;
      if (typeof newmin !== "number") newmin = 0;
      let min = Infinity,
        max = -Infinity;
      this.each(d => {
        if (d.value < min) min = d.value;
        if (d.value > max) max = d.value;
      });
      let ratio = (newmax - newmin) / (max - min);
      return this.each(d => (d.value = (d.value - min) * ratio + newmin));
    };

    /**
     * Gets the path from this Branch to `target`. If this Branch and `target` are
     * the same, returns an array containing only the Branch on which it is called.
     * @param  {Branch} target A Branch object
     * @return {Array} An ordered Array of Branches following the path between this
     * Branch and `target`
     */
    Branch.prototype.path = function(target) {
      let current = this;
      let branches = [this];
      let mrca = this.getMRCA(target);
      while (current !== mrca) {
        current = current.parent;
        branches.push(current);
      }
      let k = branches.length;
      current = target;
      while (current !== mrca) {
        branches.splice(k, 0, current);
        current = current.parent;
      }
      return branches;
    };

    /**
     * Removes a Branch and its subtree from the tree. Similar to
     * [Branch.isolate](#isolate), only it returns the root Branch of the tree
     * from which this Branch is removed.
     * @return {Branch} The root of the remaining tree.
     */
    Branch.prototype.remove = function() {
      let root = this.getRoot();
      this.isolate();
      return root;
    };

    /**
     * Removes a Branch and its subtree from the tree, and replaces it.
     * @param {Branch} replacement - The branch to replace the branch on which the
     * method is called.
     * @return {Branch} The root of the modified tree.
     */
    Branch.prototype.replace = function(replacement) {
      let root = this.getRoot();
      this.parent;
      let index = this.parent.children.indexOf(this);
      this.parent.children.splice(index, 1, replacement);
      return root;
    };

    /**
     * Reroots a tree on this Branch. Use with caution, this returns the new root,
     * which should typically supplant the existing root Branch object, but does
     * not replace that root automatically.
     * @example
     * tree = tree.children[0].children[0].reroot();
     * @return {Branch} The new root Branch, which is the Branch on which this was
     * called
     */
    Branch.prototype.reroot = function() {
      let current = this;
      let toInvert = [];
      while (!current.isRoot()) {
        toInvert.push(current);
        current = current.parent;
      }
      toInvert.reverse().forEach(c => c.invert());
      return this.fixDistances();
    };

    /**
     * Reverses the order of the children of the branch on which it is called.
     * @return {Branch} The Branch on which this was called.
     */
    Branch.prototype.rotate = function(recursive) {
      if (!this.children) return this;
      this.children.reverse();
      return this;
    };

    /**
     * Set the length of a Branch
     * @param  {number} length The new length to assign to the Branch
     * @return {Branch} The Branch object on which this was called
     */
    Branch.prototype.setLength = function(length) {
      this.length = length;
      return this;
    };

    /**
     * Sets the parent of the Branch on which it is called.
     * @param  {Branch} parent The Branch to set as parent
     * @return {Branch} The Branch on which this method was called.
     */
    Branch.prototype.setParent = function(parent) {
      if (!parent instanceof Branch && parent !== null)
        throw Error("Cannot set parent to non-Branch object!");
      this.parent = parent;
      return this;
    };

    /**
     * Collapses each descendant Branch with exactly one child into a single
     * continuous branch.
     * @return {Branch} The Branch on which this method was called.
     */
    Branch.prototype.simplify = function() {
      this.eachAfter(branch => {
        if(branch.children.length == 1){
          let child = branch.children[0];
          if(child.id == ''){
            child.id = branch.id;
          } else {
            child.id = branch.id + "+" + child.id;
          }
          branch.excise();
        }
      });
      return this.fixDistances();
    };

    /**
     * Sorts the Tree from the branch on which it is called downward.
     * @param  {Function} [comparator] A Function taking two Branches and returning
     * a numberic value. For details, see [MDN Array.sort](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/sort#Description)
     * @return {Branch} The Branch on which it was called
     */
    Branch.prototype.sort = function(comparator) {
      if (!comparator) comparator = (a, b) => a.value - b.value;
      return this.eachBefore(d => d.children.sort(comparator));
    };

    /**
     * Determines whether this Branch is likelier to be a source of `cousin`, or
     * if `cousin` is a source of this Branch.
     * @param  {Branch} cousin The other Branch to test
     * @return {Boolean} True if this might be the source of cousin, otherwise
     * false.
     */
    Branch.prototype.sources = function(cousin) {
      let mrca = this.getMRCA(cousin);
      return mrca.depthOf(this) < mrca.depthOf(cousin);
    };

    /**
     * Computes the value of each Branch according to some valuator function
     * @param  {Function} value A Function taking a Branch and returning a
     * (numeric?) value.
     * @return {Branch} The Branch on which it was called.
     */
    Branch.prototype.sum = function(value) {
      if (!value) value = d => d.value;
      return this.eachAfter(
        d => (d.value = value(d) + d.children.reduce((a, c) => a + c.value, 0))
      );
    };

    /**
     * Determines whether this Branch is likelier to be a target of `cousin`, or
     * if `cousin` is a target of this Branch.
     * @param  {Branch} cousin The other Branch to test
     * @return {Boolean} True if this might be the target of cousin, otherwise
     * false.
     */
    Branch.prototype.targets = function(cousin) {
      return cousin.sources(this);
    };

    /**
     * toJSON is an alias for [toObject](#toObject), enabling the safe use of
     * `JSON.stringify` on Branch objects (in spite of their circular references).
     * @type {Function}
     * @returns {Object} A serializable Object
     */
    Branch.prototype.toJSON = function() {
      return this.toObject();
    };

    /**
     * Computes a matrix of all patristic distances between all leaves which are
     * descendants of the Branch on which this method is called.
     * @return {Object} An Object containing a matrix (an Array of Arrays) and
     * Array of `id`s corresponding to the rows (and columns) of the matrix.
     */
    Branch.prototype.toMatrix = function() {
      let leafs = this.getLeaves();
      let n = leafs.length;
      let matrix = new Array(n);
      for (let i = 0; i < n; i++) {
        matrix[i] = new Array(n);
        matrix[i][i] = 0;
        for (let j = 0; j < i; j++) {
          let distance = leafs[i].distanceTo(leafs[j]);
          matrix[i][j] = distance;
          matrix[j][i] = distance;
        }
      }
      return {
        matrix: matrix,
        ids: leafs.map(d => d.id)
      };
    };

    /**
     * Returns the Newick representation of this Branch and its descendants.
     * @param  {Boolean} [nonterminus=falsy] Is this not the terminus of the
     * Newick Tree? This should be falsy when called by a user (i.e. you). It's
     * used internally to decide whether or not in include a semicolon in the
     * returned string.
     * @return {String} The [Newick](https://en.wikipedia.org/wiki/Newick_format)
     * representation of the Branch.
     */
    Branch.prototype.toNewick = function(nonterminus) {
      let out = "";
      if (!this.isLeaf()) {
        out +=
          "(" + this.children.map(child => child.toNewick(true)).join(",") + ")";
      }
      out += this.id;
      if (this.length) out += ":" + numberToString(this.length);
      if (!nonterminus) out += ";";
      return out;
    };

    //This function takes a number and returns a string representation that does
    //not use Scientific Notation.
    //It's adapted from [StackOverflow](https://stackoverflow.com/a/46545519/521121),
    //Which makes it available under the [CC BY-SA 3.0 License](https://creativecommons.org/licenses/by-sa/3.0/)
    function numberToString(num) {
      let numStr = String(num);
      if (Math.abs(num) < 1.0) {
        let e = parseInt(num.toString().split("e-")[1]);
        if (e) {
          let negative = num < 0;
          if (negative) num *= -1;
          num *= Math.pow(10, e - 1);
          numStr = "0." + new Array(e).join("0") + num.toString().substring(2);
          if (negative) numStr = "-" + numStr;
        }
      } else {
        let e = parseInt(num.toString().split("+")[1]);
        if (e > 20) {
          e -= 20;
          num /= Math.pow(10, e);
          numStr = num.toString() + new Array(e + 1).join("0");
        }
      }
      return numStr;
    }

    /**
     * Returns a simple Javascript object version of this Branch and its
     * descendants. This is useful in cases where you want to serialize the tree
     * (e.g. `JSON.stringify(tree)`) but can't because the tree contains circular
     * references (for simplicity, elegance, and performance reasons, each Branch
     * tracks both its children and its parent).
     * @return {Object} A serializable bare Javascript Object representing this
     * Branch and its descendants.
     */
    Branch.prototype.toObject = function() {
      let output = {
        id: this.id,
        length: this.length
      };
      if (this.children.length > 0)
        output.children = this.children.map(c => c.toObject());
      return output;
    };

    /**
     * Returns a valid JSON-string version of this Branch and its descendants.
     * @param {Function} replacer A replacer function to [pass to `JSON.stringify`](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/JSON/stringify#Parameters).
     * @param {(Number|String)} space A string or number of spaces to use for
     * indenting the output. See https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/JSON/stringify#Parameters
     * for additional details.
     * @return {Object} A valid JSON string representing this Branch and its
     * descendants.
     */
    Branch.prototype.toString = function(replacer, width) {
      if (!replacer) replacer = null;
      if (!width) width = 0;
      return JSON.stringify(this, replacer, width);
    };

    /**
     * Parses a hierarchical JSON string (or Object) as a Branch object.
     * @param  {(String|Object)} json A json string (or Javascript Object)
     * representing hierarchical data.
     * @param  {String} [idLabel="id"] The key used in the objects of `json` to
     * indicate their identifiers.
     * @param  {String} [lengthLabel='length'] The key used in the objects of `json`
     * to indicate their length.
     * @param  {String} [childrenLabel=`children`] The key used in the objects of
     * `json` to indicate their children.
     * @return {Branch} The Branch representing the root of the hierarchy
     * represented by `json`.
     */
    function parseJSON(json, idLabel, lengthLabel, childrenLabel) {
      if (!idLabel) idLabel = "id";
      if (!lengthLabel) lengthLabel = "length";
      if (!childrenLabel) childrenLabel = "children";
      if (typeof json === "string") json = JSON.parse(json);
      let root = new Branch({
        id: json[idLabel],
        length: json[lengthLabel]
      });
      if (json[childrenLabel] instanceof Array) {
        json[childrenLabel].forEach(child => {
          root.addChild(parseJSON(child));
        });
      }
      return root.fixDistances();
    }

    /**
     * Parses a matrix of distances and returns the root Branch of the output tree.
     * This is adapted from Maciej Korzepa's [neighbor-joining](https://github.com/biosustain/neighbor-joining),
     * which is released for modification under the [MIT License](https://opensource.org/licenses/MIT).
     * @param  {Array} matrix An array of `n` arrays of length `n`
     * @param  {Array} labels An array of `n` strings, each corresponding to the
     * values in `matrix`.
     * @return {Branch} A Branch object representing the root Branch of the tree
     * inferred by neighbor joining on `matrix`.
     */
    function parseMatrix(matrix, labels) {
      let that = {};
      let N = (that.N = matrix.length);
      if (!labels) labels = [...Array(N).keys()];
      that.cN = that.N;
      that.D = matrix;
      that.labels = labels;
      that.labelToTaxon = {};
      that.currIndexToLabel = new Array(N);
      that.rowChange = new Array(N);
      that.newRow = new Array(N);
      that.labelToNode = new Array(2 * N);
      that.nextIndex = N;
      that.I = new Array(that.N);
      that.S = new Array(that.N);
      for (let i = 0; i < that.N; i++) {
        let sortedRow = sortWithIndices(that.D[i], i);
        that.S[i] = sortedRow;
        that.I[i] = sortedRow.sortIndices;
      }
      that.removedIndices = new Set();
      that.indicesLeft = new Set();
      for (let i = 0; i < N; i++) {
        that.currIndexToLabel[i] = i;
        that.indicesLeft.add(i);
      }
      that.rowSumMax = 0;
      that.PNewick = "";
      let minI, minJ, d1, d2, l1, l2, node1, node2, node3;

      function setUpNode(labelIndex, distance) {
        let node;
        if (labelIndex < that.N) {
          node = new Branch({ id: that.labels[labelIndex], length: distance });
          that.labelToNode[labelIndex] = node;
        } else {
          node = that.labelToNode[labelIndex];
          node.setLength(distance);
        }
        return node;
      }

      that.rowSums = sumRows(that.D);
      for (let i = 0; i < that.cN; i++) {
        if (that.rowSums[i] > that.rowSumMax) that.rowSumMax = that.rowSums[i];
      }

      while (that.cN > 2) {
        //if (that.cN % 100 == 0 ) console.log(that.cN);
        ({ minI, minJ } = search(that));

        d1 =
          0.5 * that.D[minI][minJ] +
          (that.rowSums[minI] - that.rowSums[minJ]) / (2 * that.cN - 4);
        d2 = that.D[minI][minJ] - d1;

        l1 = that.currIndexToLabel[minI];
        l2 = that.currIndexToLabel[minJ];

        node1 = setUpNode(l1, d1);
        node2 = setUpNode(l2, d2);
        node3 = new Branch({ children: [node1, node2] });

        recalculateDistanceMatrix(that, minI, minJ);
        let sorted = sortWithIndices(that.D[minJ], minJ);
        that.S[minJ] = sorted;
        that.I[minJ] = sorted.sortIndices;
        that.S[minI] = that.I[minI] = [];
        that.cN--;

        that.labelToNode[that.nextIndex] = node3;
        that.currIndexToLabel[minI] = -1;
        that.currIndexToLabel[minJ] = that.nextIndex++;
      }

      let left = that.indicesLeft.values();
      minI = left.next().value;
      minJ = left.next().value;

      l1 = that.currIndexToLabel[minI];
      l2 = that.currIndexToLabel[minJ];
      d1 = d2 = that.D[minI][minJ] / 2;

      node1 = setUpNode(l1, d1);
      node2 = setUpNode(l2, d2);

      let tree = new Branch({ children: [node1, node2] });
      tree.fixParenthood();
      return tree.fixDistances();
    }

    function search(t) {
      let qMin = Infinity,
        D = t.D,
        cN = t.cN,
        n2 = cN - 2,
        S = t.S,
        I = t.I,
        rowSums = t.rowSums,
        removedColumns = t.removedIndices,
        uMax = t.rowSumMax,
        q,
        minI = -1,
        minJ = -1,
        c2;

      // initial guess for qMin
      for (let r = 0; r < t.N; r++) {
        if (removedColumns.has(r)) continue;
        c2 = I[r][0];
        if (removedColumns.has(c2)) continue;
        q = D[r][c2] * n2 - rowSums[r] - rowSums[c2];
        if (q < qMin) {
          qMin = q;
          minI = r;
          minJ = c2;
        }
      }

      for (let r = 0; r < t.N; r++) {
        if (removedColumns.has(r)) continue;
        for (let c = 0; c < S[r].length; c++) {
          c2 = I[r][c];
          if (removedColumns.has(c2)) continue;
          if (S[r][c] * n2 - rowSums[r] - uMax > qMin) break;
          q = D[r][c2] * n2 - rowSums[r] - rowSums[c2];
          if (q < qMin) {
            qMin = q;
            minI = r;
            minJ = c2;
          }
        }
      }

      return { minI, minJ };
    }

    function recalculateDistanceMatrix(t, joinedIndex1, joinedIndex2) {
      let D = t.D,
        n = D.length,
        sum = 0,
        aux,
        aux2,
        removedIndices = t.removedIndices,
        rowSums = t.rowSums,
        newRow = t.newRow,
        rowChange = t.rowChange,
        newMax = 0;

      removedIndices.add(joinedIndex1);
      for (let i = 0; i < n; i++) {
        if (removedIndices.has(i)) continue;
        aux = D[joinedIndex1][i] + D[joinedIndex2][i];
        aux2 = D[joinedIndex1][joinedIndex2];
        newRow[i] = 0.5 * (aux - aux2);
        sum += newRow[i];
        rowChange[i] = -0.5 * (aux + aux2);
      }
      for (let i = 0; i < n; i++) {
        D[joinedIndex1][i] = -1;
        D[i][joinedIndex1] = -1;
        if (removedIndices.has(i)) continue;
        D[joinedIndex2][i] = newRow[i];
        D[i][joinedIndex2] = newRow[i];
        rowSums[i] += rowChange[i];
        if (rowSums[i] > newMax) newMax = rowSums[i];
      }
      rowSums[joinedIndex1] = 0;
      rowSums[joinedIndex2] = sum;
      if (sum > newMax) newMax = sum;
      t.rowSumMax = newMax;
      t.indicesLeft.delete(joinedIndex1);
    }

    function sumRows(a) {
      let n = a.length,
        sums = new Array(n);
      for (let i = 0; i < n; i++) {
        let sum = 0;
        for (let j = 0; j < n; j++) {
          let v = parseFloat(a[i][j]);
          if (typeof v !== "number") continue;
          sum += a[i][j];
        }
        sums[i] = sum;
      }
      return sums;
    }

    function sortWithIndices(toSort, skip) {
      if (typeof skip === "undefined") skip = -1;
      let n = toSort.length;
      let indexCopy = new Array(n);
      let valueCopy = new Array(n);
      let i2 = 0;
      for (let i = 0; i < n; i++) {
        if (toSort[i] === -1 || i === skip) continue;
        indexCopy[i2] = i;
        valueCopy[i2++] = toSort[i];
      }
      indexCopy.length = i2;
      valueCopy.length = i2;
      indexCopy.sort((a, b) => toSort[a] - toSort[b]);
      valueCopy.sortIndices = indexCopy;
      for (let j = 0; j < i2; j++) {
        valueCopy[j] = toSort[indexCopy[j]];
      }
      return valueCopy;
    }

    /**
     * Parses a Newick String and returns a Branch object representing the root
     * of the output Tree.
     * This is adapted Jason Davies' [newick.js](https://github.com/jasondavies/newick.js/blob/master/src/newick.js),
     * which is released for modification under [the MIT License](https://opensource.org/licenses/MIT).
     * @param  {string} newick A Newick String
     * @return {Branch} A Branch representing the root of the output tree
     */
    function parseNewick(newick) {
      let ancestors = [],
        tree = new Branch(),
        tokens = newick.split(/\s*(;|\(|\)|,|:)\s*/),
        n = tokens.length;
      for (let t = 0; t < n; t++) {
        let token = tokens[t];
        let c;
        switch (token) {
          case "(": // new Branchset
            c = tree.addChild();
            ancestors.push(tree);
            tree = c;
            break;
          case ",": // another Branch
            c = ancestors[ancestors.length - 1].addChild();
            tree = c;
            break;
          case ")": // optional name next
            tree = ancestors.pop();
            break;
          case ":": // optional length next
            break;
          default:
            let x = tokens[t - 1];
            if (x == ")" || x == "(" || x == ",") {
              tree.id = token;
            } else if (x == ":") {
              tree.length = parseFloat(token);
            }
        }
      }
      return tree.fixDistances();
    }

    exports.Branch = Branch;
    exports.parseJSON = parseJSON;
    exports.parseMatrix = parseMatrix;
    exports.parseNewick = parseNewick;
    exports.version = version;

    Object.defineProperty(exports, '__esModule', { value: true });

  })));

  function identity(x) {
    return x;
  }

  var top = 1,
      right = 2,
      bottom = 3,
      left = 4,
      epsilon = 1e-6;

  function translateX(x) {
    return "translate(" + x + ",0)";
  }

  function translateY(y) {
    return "translate(0," + y + ")";
  }

  function number(scale) {
    return d => +scale(d);
  }

  function center(scale, offset) {
    offset = Math.max(0, scale.bandwidth() - offset * 2) / 2;
    if (scale.round()) offset = Math.round(offset);
    return d => +scale(d) + offset;
  }

  function entering() {
    return !this.__axis;
  }

  function axis(orient, scale) {
    var tickArguments = [],
        tickValues = null,
        tickFormat = null,
        tickSizeInner = 6,
        tickSizeOuter = 6,
        tickPadding = 3,
        offset = typeof window !== "undefined" && window.devicePixelRatio > 1 ? 0 : 0.5,
        k = orient === top || orient === left ? -1 : 1,
        x = orient === left || orient === right ? "x" : "y",
        transform = orient === top || orient === bottom ? translateX : translateY;

    function axis(context) {
      var values = tickValues == null ? (scale.ticks ? scale.ticks.apply(scale, tickArguments) : scale.domain()) : tickValues,
          format = tickFormat == null ? (scale.tickFormat ? scale.tickFormat.apply(scale, tickArguments) : identity) : tickFormat,
          spacing = Math.max(tickSizeInner, 0) + tickPadding,
          range = scale.range(),
          range0 = +range[0] + offset,
          range1 = +range[range.length - 1] + offset,
          position = (scale.bandwidth ? center : number)(scale.copy(), offset),
          selection = context.selection ? context.selection() : context,
          path = selection.selectAll(".domain").data([null]),
          tick = selection.selectAll(".tick").data(values, scale).order(),
          tickExit = tick.exit(),
          tickEnter = tick.enter().append("g").attr("class", "tick"),
          line = tick.select("line"),
          text = tick.select("text");

      path = path.merge(path.enter().insert("path", ".tick")
          .attr("class", "domain")
          .attr("stroke", "currentColor"));

      tick = tick.merge(tickEnter);

      line = line.merge(tickEnter.append("line")
          .attr("stroke", "currentColor")
          .attr(x + "2", k * tickSizeInner));

      text = text.merge(tickEnter.append("text")
          .attr("fill", "currentColor")
          .attr(x, k * spacing)
          .attr("dy", orient === top ? "0em" : orient === bottom ? "0.71em" : "0.32em"));

      if (context !== selection) {
        path = path.transition(context);
        tick = tick.transition(context);
        line = line.transition(context);
        text = text.transition(context);

        tickExit = tickExit.transition(context)
            .attr("opacity", epsilon)
            .attr("transform", function(d) { return isFinite(d = position(d)) ? transform(d + offset) : this.getAttribute("transform"); });

        tickEnter
            .attr("opacity", epsilon)
            .attr("transform", function(d) { var p = this.parentNode.__axis; return transform((p && isFinite(p = p(d)) ? p : position(d)) + offset); });
      }

      tickExit.remove();

      path
          .attr("d", orient === left || orient === right
              ? (tickSizeOuter ? "M" + k * tickSizeOuter + "," + range0 + "H" + offset + "V" + range1 + "H" + k * tickSizeOuter : "M" + offset + "," + range0 + "V" + range1)
              : (tickSizeOuter ? "M" + range0 + "," + k * tickSizeOuter + "V" + offset + "H" + range1 + "V" + k * tickSizeOuter : "M" + range0 + "," + offset + "H" + range1));

      tick
          .attr("opacity", 1)
          .attr("transform", function(d) { return transform(position(d) + offset); });

      line
          .attr(x + "2", k * tickSizeInner);

      text
          .attr(x, k * spacing)
          .text(format);

      selection.filter(entering)
          .attr("fill", "none")
          .attr("font-size", 10)
          .attr("font-family", "sans-serif")
          .attr("text-anchor", orient === right ? "start" : orient === left ? "end" : "middle");

      selection
          .each(function() { this.__axis = position; });
    }

    axis.scale = function(_) {
      return arguments.length ? (scale = _, axis) : scale;
    };

    axis.ticks = function() {
      return tickArguments = Array.from(arguments), axis;
    };

    axis.tickArguments = function(_) {
      return arguments.length ? (tickArguments = _ == null ? [] : Array.from(_), axis) : tickArguments.slice();
    };

    axis.tickValues = function(_) {
      return arguments.length ? (tickValues = _ == null ? null : Array.from(_), axis) : tickValues && tickValues.slice();
    };

    axis.tickFormat = function(_) {
      return arguments.length ? (tickFormat = _, axis) : tickFormat;
    };

    axis.tickSize = function(_) {
      return arguments.length ? (tickSizeInner = tickSizeOuter = +_, axis) : tickSizeInner;
    };

    axis.tickSizeInner = function(_) {
      return arguments.length ? (tickSizeInner = +_, axis) : tickSizeInner;
    };

    axis.tickSizeOuter = function(_) {
      return arguments.length ? (tickSizeOuter = +_, axis) : tickSizeOuter;
    };

    axis.tickPadding = function(_) {
      return arguments.length ? (tickPadding = +_, axis) : tickPadding;
    };

    axis.offset = function(_) {
      return arguments.length ? (offset = +_, axis) : offset;
    };

    return axis;
  }

  function axisBottom(scale) {
    return axis(bottom, scale);
  }

  function axisLeft(scale) {
    return axis(left, scale);
  }

  function ascending(a, b) {
    return a == null || b == null ? NaN : a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
  }

  function descending(a, b) {
    return a == null || b == null ? NaN
      : b < a ? -1
      : b > a ? 1
      : b >= a ? 0
      : NaN;
  }

  function bisector(f) {
    let compare1, compare2, delta;

    // If an accessor is specified, promote it to a comparator. In this case we
    // can test whether the search value is (self-) comparable. We can’t do this
    // for a comparator (except for specific, known comparators) because we can’t
    // tell if the comparator is symmetric, and an asymmetric comparator can’t be
    // used to test whether a single value is comparable.
    if (f.length !== 2) {
      compare1 = ascending;
      compare2 = (d, x) => ascending(f(d), x);
      delta = (d, x) => f(d) - x;
    } else {
      compare1 = f === ascending || f === descending ? f : zero;
      compare2 = f;
      delta = f;
    }

    function left(a, x, lo = 0, hi = a.length) {
      if (lo < hi) {
        if (compare1(x, x) !== 0) return hi;
        do {
          const mid = (lo + hi) >>> 1;
          if (compare2(a[mid], x) < 0) lo = mid + 1;
          else hi = mid;
        } while (lo < hi);
      }
      return lo;
    }

    function right(a, x, lo = 0, hi = a.length) {
      if (lo < hi) {
        if (compare1(x, x) !== 0) return hi;
        do {
          const mid = (lo + hi) >>> 1;
          if (compare2(a[mid], x) <= 0) lo = mid + 1;
          else hi = mid;
        } while (lo < hi);
      }
      return lo;
    }

    function center(a, x, lo = 0, hi = a.length) {
      const i = left(a, x, lo, hi - 1);
      return i > lo && delta(a[i - 1], x) > -delta(a[i], x) ? i - 1 : i;
    }

    return {left, center, right};
  }

  function zero() {
    return 0;
  }

  function number$1(x) {
    return x === null ? NaN : +x;
  }

  const ascendingBisect = bisector(ascending);
  const bisectRight = ascendingBisect.right;
  const bisectCenter = bisector(number$1).center;

  var e10 = Math.sqrt(50),
      e5 = Math.sqrt(10),
      e2 = Math.sqrt(2);

  function ticks(start, stop, count) {
    var reverse,
        i = -1,
        n,
        ticks,
        step;

    stop = +stop, start = +start, count = +count;
    if (start === stop && count > 0) return [start];
    if (reverse = stop < start) n = start, start = stop, stop = n;
    if ((step = tickIncrement(start, stop, count)) === 0 || !isFinite(step)) return [];

    if (step > 0) {
      let r0 = Math.round(start / step), r1 = Math.round(stop / step);
      if (r0 * step < start) ++r0;
      if (r1 * step > stop) --r1;
      ticks = new Array(n = r1 - r0 + 1);
      while (++i < n) ticks[i] = (r0 + i) * step;
    } else {
      step = -step;
      let r0 = Math.round(start * step), r1 = Math.round(stop * step);
      if (r0 / step < start) ++r0;
      if (r1 / step > stop) --r1;
      ticks = new Array(n = r1 - r0 + 1);
      while (++i < n) ticks[i] = (r0 + i) / step;
    }

    if (reverse) ticks.reverse();

    return ticks;
  }

  function tickIncrement(start, stop, count) {
    var step = (stop - start) / Math.max(0, count),
        power = Math.floor(Math.log(step) / Math.LN10),
        error = step / Math.pow(10, power);
    return power >= 0
        ? (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1) * Math.pow(10, power)
        : -Math.pow(10, -power) / (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1);
  }

  function tickStep(start, stop, count) {
    var step0 = Math.abs(stop - start) / Math.max(0, count),
        step1 = Math.pow(10, Math.floor(Math.log(step0) / Math.LN10)),
        error = step0 / step1;
    if (error >= e10) step1 *= 10;
    else if (error >= e5) step1 *= 5;
    else if (error >= e2) step1 *= 2;
    return stop < start ? -step1 : step1;
  }

  function initRange(domain, range) {
    switch (arguments.length) {
      case 0: break;
      case 1: this.range(domain); break;
      default: this.range(range).domain(domain); break;
    }
    return this;
  }

  function define$1(constructor, factory, prototype) {
    constructor.prototype = factory.prototype = prototype;
    prototype.constructor = constructor;
  }

  function extend(parent, definition) {
    var prototype = Object.create(parent.prototype);
    for (var key in definition) prototype[key] = definition[key];
    return prototype;
  }

  function Color() {}

  var darker = 0.7;
  var brighter = 1 / darker;

  var reI = "\\s*([+-]?\\d+)\\s*",
      reN = "\\s*([+-]?\\d*\\.?\\d+(?:[eE][+-]?\\d+)?)\\s*",
      reP = "\\s*([+-]?\\d*\\.?\\d+(?:[eE][+-]?\\d+)?)%\\s*",
      reHex = /^#([0-9a-f]{3,8})$/,
      reRgbInteger = new RegExp("^rgb\\(" + [reI, reI, reI] + "\\)$"),
      reRgbPercent = new RegExp("^rgb\\(" + [reP, reP, reP] + "\\)$"),
      reRgbaInteger = new RegExp("^rgba\\(" + [reI, reI, reI, reN] + "\\)$"),
      reRgbaPercent = new RegExp("^rgba\\(" + [reP, reP, reP, reN] + "\\)$"),
      reHslPercent = new RegExp("^hsl\\(" + [reN, reP, reP] + "\\)$"),
      reHslaPercent = new RegExp("^hsla\\(" + [reN, reP, reP, reN] + "\\)$");

  var named = {
    aliceblue: 0xf0f8ff,
    antiquewhite: 0xfaebd7,
    aqua: 0x00ffff,
    aquamarine: 0x7fffd4,
    azure: 0xf0ffff,
    beige: 0xf5f5dc,
    bisque: 0xffe4c4,
    black: 0x000000,
    blanchedalmond: 0xffebcd,
    blue: 0x0000ff,
    blueviolet: 0x8a2be2,
    brown: 0xa52a2a,
    burlywood: 0xdeb887,
    cadetblue: 0x5f9ea0,
    chartreuse: 0x7fff00,
    chocolate: 0xd2691e,
    coral: 0xff7f50,
    cornflowerblue: 0x6495ed,
    cornsilk: 0xfff8dc,
    crimson: 0xdc143c,
    cyan: 0x00ffff,
    darkblue: 0x00008b,
    darkcyan: 0x008b8b,
    darkgoldenrod: 0xb8860b,
    darkgray: 0xa9a9a9,
    darkgreen: 0x006400,
    darkgrey: 0xa9a9a9,
    darkkhaki: 0xbdb76b,
    darkmagenta: 0x8b008b,
    darkolivegreen: 0x556b2f,
    darkorange: 0xff8c00,
    darkorchid: 0x9932cc,
    darkred: 0x8b0000,
    darksalmon: 0xe9967a,
    darkseagreen: 0x8fbc8f,
    darkslateblue: 0x483d8b,
    darkslategray: 0x2f4f4f,
    darkslategrey: 0x2f4f4f,
    darkturquoise: 0x00ced1,
    darkviolet: 0x9400d3,
    deeppink: 0xff1493,
    deepskyblue: 0x00bfff,
    dimgray: 0x696969,
    dimgrey: 0x696969,
    dodgerblue: 0x1e90ff,
    firebrick: 0xb22222,
    floralwhite: 0xfffaf0,
    forestgreen: 0x228b22,
    fuchsia: 0xff00ff,
    gainsboro: 0xdcdcdc,
    ghostwhite: 0xf8f8ff,
    gold: 0xffd700,
    goldenrod: 0xdaa520,
    gray: 0x808080,
    green: 0x008000,
    greenyellow: 0xadff2f,
    grey: 0x808080,
    honeydew: 0xf0fff0,
    hotpink: 0xff69b4,
    indianred: 0xcd5c5c,
    indigo: 0x4b0082,
    ivory: 0xfffff0,
    khaki: 0xf0e68c,
    lavender: 0xe6e6fa,
    lavenderblush: 0xfff0f5,
    lawngreen: 0x7cfc00,
    lemonchiffon: 0xfffacd,
    lightblue: 0xadd8e6,
    lightcoral: 0xf08080,
    lightcyan: 0xe0ffff,
    lightgoldenrodyellow: 0xfafad2,
    lightgray: 0xd3d3d3,
    lightgreen: 0x90ee90,
    lightgrey: 0xd3d3d3,
    lightpink: 0xffb6c1,
    lightsalmon: 0xffa07a,
    lightseagreen: 0x20b2aa,
    lightskyblue: 0x87cefa,
    lightslategray: 0x778899,
    lightslategrey: 0x778899,
    lightsteelblue: 0xb0c4de,
    lightyellow: 0xffffe0,
    lime: 0x00ff00,
    limegreen: 0x32cd32,
    linen: 0xfaf0e6,
    magenta: 0xff00ff,
    maroon: 0x800000,
    mediumaquamarine: 0x66cdaa,
    mediumblue: 0x0000cd,
    mediumorchid: 0xba55d3,
    mediumpurple: 0x9370db,
    mediumseagreen: 0x3cb371,
    mediumslateblue: 0x7b68ee,
    mediumspringgreen: 0x00fa9a,
    mediumturquoise: 0x48d1cc,
    mediumvioletred: 0xc71585,
    midnightblue: 0x191970,
    mintcream: 0xf5fffa,
    mistyrose: 0xffe4e1,
    moccasin: 0xffe4b5,
    navajowhite: 0xffdead,
    navy: 0x000080,
    oldlace: 0xfdf5e6,
    olive: 0x808000,
    olivedrab: 0x6b8e23,
    orange: 0xffa500,
    orangered: 0xff4500,
    orchid: 0xda70d6,
    palegoldenrod: 0xeee8aa,
    palegreen: 0x98fb98,
    paleturquoise: 0xafeeee,
    palevioletred: 0xdb7093,
    papayawhip: 0xffefd5,
    peachpuff: 0xffdab9,
    peru: 0xcd853f,
    pink: 0xffc0cb,
    plum: 0xdda0dd,
    powderblue: 0xb0e0e6,
    purple: 0x800080,
    rebeccapurple: 0x663399,
    red: 0xff0000,
    rosybrown: 0xbc8f8f,
    royalblue: 0x4169e1,
    saddlebrown: 0x8b4513,
    salmon: 0xfa8072,
    sandybrown: 0xf4a460,
    seagreen: 0x2e8b57,
    seashell: 0xfff5ee,
    sienna: 0xa0522d,
    silver: 0xc0c0c0,
    skyblue: 0x87ceeb,
    slateblue: 0x6a5acd,
    slategray: 0x708090,
    slategrey: 0x708090,
    snow: 0xfffafa,
    springgreen: 0x00ff7f,
    steelblue: 0x4682b4,
    tan: 0xd2b48c,
    teal: 0x008080,
    thistle: 0xd8bfd8,
    tomato: 0xff6347,
    turquoise: 0x40e0d0,
    violet: 0xee82ee,
    wheat: 0xf5deb3,
    white: 0xffffff,
    whitesmoke: 0xf5f5f5,
    yellow: 0xffff00,
    yellowgreen: 0x9acd32
  };

  define$1(Color, color, {
    copy: function(channels) {
      return Object.assign(new this.constructor, this, channels);
    },
    displayable: function() {
      return this.rgb().displayable();
    },
    hex: color_formatHex, // Deprecated! Use color.formatHex.
    formatHex: color_formatHex,
    formatHsl: color_formatHsl,
    formatRgb: color_formatRgb,
    toString: color_formatRgb
  });

  function color_formatHex() {
    return this.rgb().formatHex();
  }

  function color_formatHsl() {
    return hslConvert(this).formatHsl();
  }

  function color_formatRgb() {
    return this.rgb().formatRgb();
  }

  function color(format) {
    var m, l;
    format = (format + "").trim().toLowerCase();
    return (m = reHex.exec(format)) ? (l = m[1].length, m = parseInt(m[1], 16), l === 6 ? rgbn(m) // #ff0000
        : l === 3 ? new Rgb((m >> 8 & 0xf) | (m >> 4 & 0xf0), (m >> 4 & 0xf) | (m & 0xf0), ((m & 0xf) << 4) | (m & 0xf), 1) // #f00
        : l === 8 ? rgba(m >> 24 & 0xff, m >> 16 & 0xff, m >> 8 & 0xff, (m & 0xff) / 0xff) // #ff000000
        : l === 4 ? rgba((m >> 12 & 0xf) | (m >> 8 & 0xf0), (m >> 8 & 0xf) | (m >> 4 & 0xf0), (m >> 4 & 0xf) | (m & 0xf0), (((m & 0xf) << 4) | (m & 0xf)) / 0xff) // #f000
        : null) // invalid hex
        : (m = reRgbInteger.exec(format)) ? new Rgb(m[1], m[2], m[3], 1) // rgb(255, 0, 0)
        : (m = reRgbPercent.exec(format)) ? new Rgb(m[1] * 255 / 100, m[2] * 255 / 100, m[3] * 255 / 100, 1) // rgb(100%, 0%, 0%)
        : (m = reRgbaInteger.exec(format)) ? rgba(m[1], m[2], m[3], m[4]) // rgba(255, 0, 0, 1)
        : (m = reRgbaPercent.exec(format)) ? rgba(m[1] * 255 / 100, m[2] * 255 / 100, m[3] * 255 / 100, m[4]) // rgb(100%, 0%, 0%, 1)
        : (m = reHslPercent.exec(format)) ? hsla(m[1], m[2] / 100, m[3] / 100, 1) // hsl(120, 50%, 50%)
        : (m = reHslaPercent.exec(format)) ? hsla(m[1], m[2] / 100, m[3] / 100, m[4]) // hsla(120, 50%, 50%, 1)
        : named.hasOwnProperty(format) ? rgbn(named[format]) // eslint-disable-line no-prototype-builtins
        : format === "transparent" ? new Rgb(NaN, NaN, NaN, 0)
        : null;
  }

  function rgbn(n) {
    return new Rgb(n >> 16 & 0xff, n >> 8 & 0xff, n & 0xff, 1);
  }

  function rgba(r, g, b, a) {
    if (a <= 0) r = g = b = NaN;
    return new Rgb(r, g, b, a);
  }

  function rgbConvert(o) {
    if (!(o instanceof Color)) o = color(o);
    if (!o) return new Rgb;
    o = o.rgb();
    return new Rgb(o.r, o.g, o.b, o.opacity);
  }

  function rgb(r, g, b, opacity) {
    return arguments.length === 1 ? rgbConvert(r) : new Rgb(r, g, b, opacity == null ? 1 : opacity);
  }

  function Rgb(r, g, b, opacity) {
    this.r = +r;
    this.g = +g;
    this.b = +b;
    this.opacity = +opacity;
  }

  define$1(Rgb, rgb, extend(Color, {
    brighter: function(k) {
      k = k == null ? brighter : Math.pow(brighter, k);
      return new Rgb(this.r * k, this.g * k, this.b * k, this.opacity);
    },
    darker: function(k) {
      k = k == null ? darker : Math.pow(darker, k);
      return new Rgb(this.r * k, this.g * k, this.b * k, this.opacity);
    },
    rgb: function() {
      return this;
    },
    displayable: function() {
      return (-0.5 <= this.r && this.r < 255.5)
          && (-0.5 <= this.g && this.g < 255.5)
          && (-0.5 <= this.b && this.b < 255.5)
          && (0 <= this.opacity && this.opacity <= 1);
    },
    hex: rgb_formatHex, // Deprecated! Use color.formatHex.
    formatHex: rgb_formatHex,
    formatRgb: rgb_formatRgb,
    toString: rgb_formatRgb
  }));

  function rgb_formatHex() {
    return "#" + hex(this.r) + hex(this.g) + hex(this.b);
  }

  function rgb_formatRgb() {
    var a = this.opacity; a = isNaN(a) ? 1 : Math.max(0, Math.min(1, a));
    return (a === 1 ? "rgb(" : "rgba(")
        + Math.max(0, Math.min(255, Math.round(this.r) || 0)) + ", "
        + Math.max(0, Math.min(255, Math.round(this.g) || 0)) + ", "
        + Math.max(0, Math.min(255, Math.round(this.b) || 0))
        + (a === 1 ? ")" : ", " + a + ")");
  }

  function hex(value) {
    value = Math.max(0, Math.min(255, Math.round(value) || 0));
    return (value < 16 ? "0" : "") + value.toString(16);
  }

  function hsla(h, s, l, a) {
    if (a <= 0) h = s = l = NaN;
    else if (l <= 0 || l >= 1) h = s = NaN;
    else if (s <= 0) h = NaN;
    return new Hsl(h, s, l, a);
  }

  function hslConvert(o) {
    if (o instanceof Hsl) return new Hsl(o.h, o.s, o.l, o.opacity);
    if (!(o instanceof Color)) o = color(o);
    if (!o) return new Hsl;
    if (o instanceof Hsl) return o;
    o = o.rgb();
    var r = o.r / 255,
        g = o.g / 255,
        b = o.b / 255,
        min = Math.min(r, g, b),
        max = Math.max(r, g, b),
        h = NaN,
        s = max - min,
        l = (max + min) / 2;
    if (s) {
      if (r === max) h = (g - b) / s + (g < b) * 6;
      else if (g === max) h = (b - r) / s + 2;
      else h = (r - g) / s + 4;
      s /= l < 0.5 ? max + min : 2 - max - min;
      h *= 60;
    } else {
      s = l > 0 && l < 1 ? 0 : h;
    }
    return new Hsl(h, s, l, o.opacity);
  }

  function hsl(h, s, l, opacity) {
    return arguments.length === 1 ? hslConvert(h) : new Hsl(h, s, l, opacity == null ? 1 : opacity);
  }

  function Hsl(h, s, l, opacity) {
    this.h = +h;
    this.s = +s;
    this.l = +l;
    this.opacity = +opacity;
  }

  define$1(Hsl, hsl, extend(Color, {
    brighter: function(k) {
      k = k == null ? brighter : Math.pow(brighter, k);
      return new Hsl(this.h, this.s, this.l * k, this.opacity);
    },
    darker: function(k) {
      k = k == null ? darker : Math.pow(darker, k);
      return new Hsl(this.h, this.s, this.l * k, this.opacity);
    },
    rgb: function() {
      var h = this.h % 360 + (this.h < 0) * 360,
          s = isNaN(h) || isNaN(this.s) ? 0 : this.s,
          l = this.l,
          m2 = l + (l < 0.5 ? l : 1 - l) * s,
          m1 = 2 * l - m2;
      return new Rgb(
        hsl2rgb(h >= 240 ? h - 240 : h + 120, m1, m2),
        hsl2rgb(h, m1, m2),
        hsl2rgb(h < 120 ? h + 240 : h - 120, m1, m2),
        this.opacity
      );
    },
    displayable: function() {
      return (0 <= this.s && this.s <= 1 || isNaN(this.s))
          && (0 <= this.l && this.l <= 1)
          && (0 <= this.opacity && this.opacity <= 1);
    },
    formatHsl: function() {
      var a = this.opacity; a = isNaN(a) ? 1 : Math.max(0, Math.min(1, a));
      return (a === 1 ? "hsl(" : "hsla(")
          + (this.h || 0) + ", "
          + (this.s || 0) * 100 + "%, "
          + (this.l || 0) * 100 + "%"
          + (a === 1 ? ")" : ", " + a + ")");
    }
  }));

  /* From FvD 13.37, CSS Color Module Level 3 */
  function hsl2rgb(h, m1, m2) {
    return (h < 60 ? m1 + (m2 - m1) * h / 60
        : h < 180 ? m2
        : h < 240 ? m1 + (m2 - m1) * (240 - h) / 60
        : m1) * 255;
  }

  function constant(x) {
    return function() {
      return x;
    };
  }

  function linear(a, d) {
    return function(t) {
      return a + t * d;
    };
  }

  function exponential(a, b, y) {
    return a = Math.pow(a, y), b = Math.pow(b, y) - a, y = 1 / y, function(t) {
      return Math.pow(a + t * b, y);
    };
  }

  function gamma(y) {
    return (y = +y) === 1 ? nogamma : function(a, b) {
      return b - a ? exponential(a, b, y) : constant(isNaN(a) ? b : a);
    };
  }

  function nogamma(a, b) {
    var d = b - a;
    return d ? linear(a, d) : constant(isNaN(a) ? b : a);
  }

  var interpolateRgb = (function rgbGamma(y) {
    var color = gamma(y);

    function rgb$1(start, end) {
      var r = color((start = rgb(start)).r, (end = rgb(end)).r),
          g = color(start.g, end.g),
          b = color(start.b, end.b),
          opacity = nogamma(start.opacity, end.opacity);
      return function(t) {
        start.r = r(t);
        start.g = g(t);
        start.b = b(t);
        start.opacity = opacity(t);
        return start + "";
      };
    }

    rgb$1.gamma = rgbGamma;

    return rgb$1;
  })(1);

  function numberArray(a, b) {
    if (!b) b = [];
    var n = a ? Math.min(b.length, a.length) : 0,
        c = b.slice(),
        i;
    return function(t) {
      for (i = 0; i < n; ++i) c[i] = a[i] * (1 - t) + b[i] * t;
      return c;
    };
  }

  function isNumberArray(x) {
    return ArrayBuffer.isView(x) && !(x instanceof DataView);
  }

  function genericArray(a, b) {
    var nb = b ? b.length : 0,
        na = a ? Math.min(nb, a.length) : 0,
        x = new Array(na),
        c = new Array(nb),
        i;

    for (i = 0; i < na; ++i) x[i] = interpolate(a[i], b[i]);
    for (; i < nb; ++i) c[i] = b[i];

    return function(t) {
      for (i = 0; i < na; ++i) c[i] = x[i](t);
      return c;
    };
  }

  function date(a, b) {
    var d = new Date;
    return a = +a, b = +b, function(t) {
      return d.setTime(a * (1 - t) + b * t), d;
    };
  }

  function interpolateNumber(a, b) {
    return a = +a, b = +b, function(t) {
      return a * (1 - t) + b * t;
    };
  }

  function object(a, b) {
    var i = {},
        c = {},
        k;

    if (a === null || typeof a !== "object") a = {};
    if (b === null || typeof b !== "object") b = {};

    for (k in b) {
      if (k in a) {
        i[k] = interpolate(a[k], b[k]);
      } else {
        c[k] = b[k];
      }
    }

    return function(t) {
      for (k in i) c[k] = i[k](t);
      return c;
    };
  }

  var reA = /[-+]?(?:\d+\.?\d*|\.?\d+)(?:[eE][-+]?\d+)?/g,
      reB = new RegExp(reA.source, "g");

  function zero$1(b) {
    return function() {
      return b;
    };
  }

  function one(b) {
    return function(t) {
      return b(t) + "";
    };
  }

  function interpolateString(a, b) {
    var bi = reA.lastIndex = reB.lastIndex = 0, // scan index for next number in b
        am, // current match in a
        bm, // current match in b
        bs, // string preceding current number in b, if any
        i = -1, // index in s
        s = [], // string constants and placeholders
        q = []; // number interpolators

    // Coerce inputs to strings.
    a = a + "", b = b + "";

    // Interpolate pairs of numbers in a & b.
    while ((am = reA.exec(a))
        && (bm = reB.exec(b))) {
      if ((bs = bm.index) > bi) { // a string precedes the next number in b
        bs = b.slice(bi, bs);
        if (s[i]) s[i] += bs; // coalesce with previous string
        else s[++i] = bs;
      }
      if ((am = am[0]) === (bm = bm[0])) { // numbers in a & b match
        if (s[i]) s[i] += bm; // coalesce with previous string
        else s[++i] = bm;
      } else { // interpolate non-matching numbers
        s[++i] = null;
        q.push({i: i, x: interpolateNumber(am, bm)});
      }
      bi = reB.lastIndex;
    }

    // Add remains of b.
    if (bi < b.length) {
      bs = b.slice(bi);
      if (s[i]) s[i] += bs; // coalesce with previous string
      else s[++i] = bs;
    }

    // Special optimization for only a single match.
    // Otherwise, interpolate each of the numbers and rejoin the string.
    return s.length < 2 ? (q[0]
        ? one(q[0].x)
        : zero$1(b))
        : (b = q.length, function(t) {
            for (var i = 0, o; i < b; ++i) s[(o = q[i]).i] = o.x(t);
            return s.join("");
          });
  }

  function interpolate(a, b) {
    var t = typeof b, c;
    return b == null || t === "boolean" ? constant(b)
        : (t === "number" ? interpolateNumber
        : t === "string" ? ((c = color(b)) ? (b = c, interpolateRgb) : interpolateString)
        : b instanceof color ? interpolateRgb
        : b instanceof Date ? date
        : isNumberArray(b) ? numberArray
        : Array.isArray(b) ? genericArray
        : typeof b.valueOf !== "function" && typeof b.toString !== "function" || isNaN(b) ? object
        : interpolateNumber)(a, b);
  }

  function interpolateRound(a, b) {
    return a = +a, b = +b, function(t) {
      return Math.round(a * (1 - t) + b * t);
    };
  }

  var degrees = 180 / Math.PI;

  var identity$1 = {
    translateX: 0,
    translateY: 0,
    rotate: 0,
    skewX: 0,
    scaleX: 1,
    scaleY: 1
  };

  function decompose(a, b, c, d, e, f) {
    var scaleX, scaleY, skewX;
    if (scaleX = Math.sqrt(a * a + b * b)) a /= scaleX, b /= scaleX;
    if (skewX = a * c + b * d) c -= a * skewX, d -= b * skewX;
    if (scaleY = Math.sqrt(c * c + d * d)) c /= scaleY, d /= scaleY, skewX /= scaleY;
    if (a * d < b * c) a = -a, b = -b, skewX = -skewX, scaleX = -scaleX;
    return {
      translateX: e,
      translateY: f,
      rotate: Math.atan2(b, a) * degrees,
      skewX: Math.atan(skewX) * degrees,
      scaleX: scaleX,
      scaleY: scaleY
    };
  }

  var cssNode,
      cssRoot,
      cssView,
      svgNode;

  function parseCss(value) {
    if (value === "none") return identity$1;
    if (!cssNode) cssNode = document.createElement("DIV"), cssRoot = document.documentElement, cssView = document.defaultView;
    cssNode.style.transform = value;
    value = cssView.getComputedStyle(cssRoot.appendChild(cssNode), null).getPropertyValue("transform");
    cssRoot.removeChild(cssNode);
    value = value.slice(7, -1).split(",");
    return decompose(+value[0], +value[1], +value[2], +value[3], +value[4], +value[5]);
  }

  function parseSvg(value) {
    if (value == null) return identity$1;
    if (!svgNode) svgNode = document.createElementNS("http://www.w3.org/2000/svg", "g");
    svgNode.setAttribute("transform", value);
    if (!(value = svgNode.transform.baseVal.consolidate())) return identity$1;
    value = value.matrix;
    return decompose(value.a, value.b, value.c, value.d, value.e, value.f);
  }

  function interpolateTransform(parse, pxComma, pxParen, degParen) {

    function pop(s) {
      return s.length ? s.pop() + " " : "";
    }

    function translate(xa, ya, xb, yb, s, q) {
      if (xa !== xb || ya !== yb) {
        var i = s.push("translate(", null, pxComma, null, pxParen);
        q.push({i: i - 4, x: interpolateNumber(xa, xb)}, {i: i - 2, x: interpolateNumber(ya, yb)});
      } else if (xb || yb) {
        s.push("translate(" + xb + pxComma + yb + pxParen);
      }
    }

    function rotate(a, b, s, q) {
      if (a !== b) {
        if (a - b > 180) b += 360; else if (b - a > 180) a += 360; // shortest path
        q.push({i: s.push(pop(s) + "rotate(", null, degParen) - 2, x: interpolateNumber(a, b)});
      } else if (b) {
        s.push(pop(s) + "rotate(" + b + degParen);
      }
    }

    function skewX(a, b, s, q) {
      if (a !== b) {
        q.push({i: s.push(pop(s) + "skewX(", null, degParen) - 2, x: interpolateNumber(a, b)});
      } else if (b) {
        s.push(pop(s) + "skewX(" + b + degParen);
      }
    }

    function scale(xa, ya, xb, yb, s, q) {
      if (xa !== xb || ya !== yb) {
        var i = s.push(pop(s) + "scale(", null, ",", null, ")");
        q.push({i: i - 4, x: interpolateNumber(xa, xb)}, {i: i - 2, x: interpolateNumber(ya, yb)});
      } else if (xb !== 1 || yb !== 1) {
        s.push(pop(s) + "scale(" + xb + "," + yb + ")");
      }
    }

    return function(a, b) {
      var s = [], // string constants and placeholders
          q = []; // number interpolators
      a = parse(a), b = parse(b);
      translate(a.translateX, a.translateY, b.translateX, b.translateY, s, q);
      rotate(a.rotate, b.rotate, s, q);
      skewX(a.skewX, b.skewX, s, q);
      scale(a.scaleX, a.scaleY, b.scaleX, b.scaleY, s, q);
      a = b = null; // gc
      return function(t) {
        var i = -1, n = q.length, o;
        while (++i < n) s[(o = q[i]).i] = o.x(t);
        return s.join("");
      };
    };
  }

  var interpolateTransformCss = interpolateTransform(parseCss, "px, ", "px)", "deg)");
  var interpolateTransformSvg = interpolateTransform(parseSvg, ", ", ")", ")");

  var rho = Math.SQRT2,
      rho2 = 2,
      rho4 = 4,
      epsilon2 = 1e-12;

  function cosh(x) {
    return ((x = Math.exp(x)) + 1 / x) / 2;
  }

  function sinh(x) {
    return ((x = Math.exp(x)) - 1 / x) / 2;
  }

  function tanh(x) {
    return ((x = Math.exp(2 * x)) - 1) / (x + 1);
  }

  // p0 = [ux0, uy0, w0]
  // p1 = [ux1, uy1, w1]
  function interpolateZoom(p0, p1) {
    var ux0 = p0[0], uy0 = p0[1], w0 = p0[2],
        ux1 = p1[0], uy1 = p1[1], w1 = p1[2],
        dx = ux1 - ux0,
        dy = uy1 - uy0,
        d2 = dx * dx + dy * dy,
        i,
        S;

    // Special case for u0 ≅ u1.
    if (d2 < epsilon2) {
      S = Math.log(w1 / w0) / rho;
      i = function(t) {
        return [
          ux0 + t * dx,
          uy0 + t * dy,
          w0 * Math.exp(rho * t * S)
        ];
      };
    }

    // General case.
    else {
      var d1 = Math.sqrt(d2),
          b0 = (w1 * w1 - w0 * w0 + rho4 * d2) / (2 * w0 * rho2 * d1),
          b1 = (w1 * w1 - w0 * w0 - rho4 * d2) / (2 * w1 * rho2 * d1),
          r0 = Math.log(Math.sqrt(b0 * b0 + 1) - b0),
          r1 = Math.log(Math.sqrt(b1 * b1 + 1) - b1);
      S = (r1 - r0) / rho;
      i = function(t) {
        var s = t * S,
            coshr0 = cosh(r0),
            u = w0 / (rho2 * d1) * (coshr0 * tanh(rho * s + r0) - sinh(r0));
        return [
          ux0 + u * dx,
          uy0 + u * dy,
          w0 * coshr0 / cosh(rho * s + r0)
        ];
      };
    }

    i.duration = S * 1000;

    return i;
  }

  function constants(x) {
    return function() {
      return x;
    };
  }

  function number$2(x) {
    return +x;
  }

  var unit = [0, 1];

  function identity$2(x) {
    return x;
  }

  function normalize(a, b) {
    return (b -= (a = +a))
        ? function(x) { return (x - a) / b; }
        : constants(isNaN(b) ? NaN : 0.5);
  }

  function clamper(a, b) {
    var t;
    if (a > b) t = a, a = b, b = t;
    return function(x) { return Math.max(a, Math.min(b, x)); };
  }

  // normalize(a, b)(x) takes a domain value x in [a,b] and returns the corresponding parameter t in [0,1].
  // interpolate(a, b)(t) takes a parameter t in [0,1] and returns the corresponding range value x in [a,b].
  function bimap(domain, range, interpolate) {
    var d0 = domain[0], d1 = domain[1], r0 = range[0], r1 = range[1];
    if (d1 < d0) d0 = normalize(d1, d0), r0 = interpolate(r1, r0);
    else d0 = normalize(d0, d1), r0 = interpolate(r0, r1);
    return function(x) { return r0(d0(x)); };
  }

  function polymap(domain, range, interpolate) {
    var j = Math.min(domain.length, range.length) - 1,
        d = new Array(j),
        r = new Array(j),
        i = -1;

    // Reverse descending domains.
    if (domain[j] < domain[0]) {
      domain = domain.slice().reverse();
      range = range.slice().reverse();
    }

    while (++i < j) {
      d[i] = normalize(domain[i], domain[i + 1]);
      r[i] = interpolate(range[i], range[i + 1]);
    }

    return function(x) {
      var i = bisectRight(domain, x, 1, j) - 1;
      return r[i](d[i](x));
    };
  }

  function copy(source, target) {
    return target
        .domain(source.domain())
        .range(source.range())
        .interpolate(source.interpolate())
        .clamp(source.clamp())
        .unknown(source.unknown());
  }

  function transformer() {
    var domain = unit,
        range = unit,
        interpolate$1 = interpolate,
        transform,
        untransform,
        unknown,
        clamp = identity$2,
        piecewise,
        output,
        input;

    function rescale() {
      var n = Math.min(domain.length, range.length);
      if (clamp !== identity$2) clamp = clamper(domain[0], domain[n - 1]);
      piecewise = n > 2 ? polymap : bimap;
      output = input = null;
      return scale;
    }

    function scale(x) {
      return x == null || isNaN(x = +x) ? unknown : (output || (output = piecewise(domain.map(transform), range, interpolate$1)))(transform(clamp(x)));
    }

    scale.invert = function(y) {
      return clamp(untransform((input || (input = piecewise(range, domain.map(transform), interpolateNumber)))(y)));
    };

    scale.domain = function(_) {
      return arguments.length ? (domain = Array.from(_, number$2), rescale()) : domain.slice();
    };

    scale.range = function(_) {
      return arguments.length ? (range = Array.from(_), rescale()) : range.slice();
    };

    scale.rangeRound = function(_) {
      return range = Array.from(_), interpolate$1 = interpolateRound, rescale();
    };

    scale.clamp = function(_) {
      return arguments.length ? (clamp = _ ? true : identity$2, rescale()) : clamp !== identity$2;
    };

    scale.interpolate = function(_) {
      return arguments.length ? (interpolate$1 = _, rescale()) : interpolate$1;
    };

    scale.unknown = function(_) {
      return arguments.length ? (unknown = _, scale) : unknown;
    };

    return function(t, u) {
      transform = t, untransform = u;
      return rescale();
    };
  }

  function continuous() {
    return transformer()(identity$2, identity$2);
  }

  function formatDecimal(x) {
    return Math.abs(x = Math.round(x)) >= 1e21
        ? x.toLocaleString("en").replace(/,/g, "")
        : x.toString(10);
  }

  // Computes the decimal coefficient and exponent of the specified number x with
  // significant digits p, where x is positive and p is in [1, 21] or undefined.
  // For example, formatDecimalParts(1.23) returns ["123", 0].
  function formatDecimalParts(x, p) {
    if ((i = (x = p ? x.toExponential(p - 1) : x.toExponential()).indexOf("e")) < 0) return null; // NaN, ±Infinity
    var i, coefficient = x.slice(0, i);

    // The string returned by toExponential either has the form \d\.\d+e[-+]\d+
    // (e.g., 1.2e+3) or the form \de[-+]\d+ (e.g., 1e+3).
    return [
      coefficient.length > 1 ? coefficient[0] + coefficient.slice(2) : coefficient,
      +x.slice(i + 1)
    ];
  }

  function exponent(x) {
    return x = formatDecimalParts(Math.abs(x)), x ? x[1] : NaN;
  }

  function formatGroup(grouping, thousands) {
    return function(value, width) {
      var i = value.length,
          t = [],
          j = 0,
          g = grouping[0],
          length = 0;

      while (i > 0 && g > 0) {
        if (length + g + 1 > width) g = Math.max(1, width - length);
        t.push(value.substring(i -= g, i + g));
        if ((length += g + 1) > width) break;
        g = grouping[j = (j + 1) % grouping.length];
      }

      return t.reverse().join(thousands);
    };
  }

  function formatNumerals(numerals) {
    return function(value) {
      return value.replace(/[0-9]/g, function(i) {
        return numerals[+i];
      });
    };
  }

  // [[fill]align][sign][symbol][0][width][,][.precision][~][type]
  var re = /^(?:(.)?([<>=^]))?([+\-( ])?([$#])?(0)?(\d+)?(,)?(\.\d+)?(~)?([a-z%])?$/i;

  function formatSpecifier(specifier) {
    if (!(match = re.exec(specifier))) throw new Error("invalid format: " + specifier);
    var match;
    return new FormatSpecifier({
      fill: match[1],
      align: match[2],
      sign: match[3],
      symbol: match[4],
      zero: match[5],
      width: match[6],
      comma: match[7],
      precision: match[8] && match[8].slice(1),
      trim: match[9],
      type: match[10]
    });
  }

  formatSpecifier.prototype = FormatSpecifier.prototype; // instanceof

  function FormatSpecifier(specifier) {
    this.fill = specifier.fill === undefined ? " " : specifier.fill + "";
    this.align = specifier.align === undefined ? ">" : specifier.align + "";
    this.sign = specifier.sign === undefined ? "-" : specifier.sign + "";
    this.symbol = specifier.symbol === undefined ? "" : specifier.symbol + "";
    this.zero = !!specifier.zero;
    this.width = specifier.width === undefined ? undefined : +specifier.width;
    this.comma = !!specifier.comma;
    this.precision = specifier.precision === undefined ? undefined : +specifier.precision;
    this.trim = !!specifier.trim;
    this.type = specifier.type === undefined ? "" : specifier.type + "";
  }

  FormatSpecifier.prototype.toString = function() {
    return this.fill
        + this.align
        + this.sign
        + this.symbol
        + (this.zero ? "0" : "")
        + (this.width === undefined ? "" : Math.max(1, this.width | 0))
        + (this.comma ? "," : "")
        + (this.precision === undefined ? "" : "." + Math.max(0, this.precision | 0))
        + (this.trim ? "~" : "")
        + this.type;
  };

  // Trims insignificant zeros, e.g., replaces 1.2000k with 1.2k.
  function formatTrim(s) {
    out: for (var n = s.length, i = 1, i0 = -1, i1; i < n; ++i) {
      switch (s[i]) {
        case ".": i0 = i1 = i; break;
        case "0": if (i0 === 0) i0 = i; i1 = i; break;
        default: if (!+s[i]) break out; if (i0 > 0) i0 = 0; break;
      }
    }
    return i0 > 0 ? s.slice(0, i0) + s.slice(i1 + 1) : s;
  }

  var prefixExponent;

  function formatPrefixAuto(x, p) {
    var d = formatDecimalParts(x, p);
    if (!d) return x + "";
    var coefficient = d[0],
        exponent = d[1],
        i = exponent - (prefixExponent = Math.max(-8, Math.min(8, Math.floor(exponent / 3))) * 3) + 1,
        n = coefficient.length;
    return i === n ? coefficient
        : i > n ? coefficient + new Array(i - n + 1).join("0")
        : i > 0 ? coefficient.slice(0, i) + "." + coefficient.slice(i)
        : "0." + new Array(1 - i).join("0") + formatDecimalParts(x, Math.max(0, p + i - 1))[0]; // less than 1y!
  }

  function formatRounded(x, p) {
    var d = formatDecimalParts(x, p);
    if (!d) return x + "";
    var coefficient = d[0],
        exponent = d[1];
    return exponent < 0 ? "0." + new Array(-exponent).join("0") + coefficient
        : coefficient.length > exponent + 1 ? coefficient.slice(0, exponent + 1) + "." + coefficient.slice(exponent + 1)
        : coefficient + new Array(exponent - coefficient.length + 2).join("0");
  }

  var formatTypes = {
    "%": function(x, p) { return (x * 100).toFixed(p); },
    "b": function(x) { return Math.round(x).toString(2); },
    "c": function(x) { return x + ""; },
    "d": formatDecimal,
    "e": function(x, p) { return x.toExponential(p); },
    "f": function(x, p) { return x.toFixed(p); },
    "g": function(x, p) { return x.toPrecision(p); },
    "o": function(x) { return Math.round(x).toString(8); },
    "p": function(x, p) { return formatRounded(x * 100, p); },
    "r": formatRounded,
    "s": formatPrefixAuto,
    "X": function(x) { return Math.round(x).toString(16).toUpperCase(); },
    "x": function(x) { return Math.round(x).toString(16); }
  };

  function identity$3(x) {
    return x;
  }

  var map = Array.prototype.map,
      prefixes = ["y","z","a","f","p","n","µ","m","","k","M","G","T","P","E","Z","Y"];

  function formatLocale(locale) {
    var group = locale.grouping === undefined || locale.thousands === undefined ? identity$3 : formatGroup(map.call(locale.grouping, Number), locale.thousands + ""),
        currencyPrefix = locale.currency === undefined ? "" : locale.currency[0] + "",
        currencySuffix = locale.currency === undefined ? "" : locale.currency[1] + "",
        decimal = locale.decimal === undefined ? "." : locale.decimal + "",
        numerals = locale.numerals === undefined ? identity$3 : formatNumerals(map.call(locale.numerals, String)),
        percent = locale.percent === undefined ? "%" : locale.percent + "",
        minus = locale.minus === undefined ? "-" : locale.minus + "",
        nan = locale.nan === undefined ? "NaN" : locale.nan + "";

    function newFormat(specifier) {
      specifier = formatSpecifier(specifier);

      var fill = specifier.fill,
          align = specifier.align,
          sign = specifier.sign,
          symbol = specifier.symbol,
          zero = specifier.zero,
          width = specifier.width,
          comma = specifier.comma,
          precision = specifier.precision,
          trim = specifier.trim,
          type = specifier.type;

      // The "n" type is an alias for ",g".
      if (type === "n") comma = true, type = "g";

      // The "" type, and any invalid type, is an alias for ".12~g".
      else if (!formatTypes[type]) precision === undefined && (precision = 12), trim = true, type = "g";

      // If zero fill is specified, padding goes after sign and before digits.
      if (zero || (fill === "0" && align === "=")) zero = true, fill = "0", align = "=";

      // Compute the prefix and suffix.
      // For SI-prefix, the suffix is lazily computed.
      var prefix = symbol === "$" ? currencyPrefix : symbol === "#" && /[boxX]/.test(type) ? "0" + type.toLowerCase() : "",
          suffix = symbol === "$" ? currencySuffix : /[%p]/.test(type) ? percent : "";

      // What format function should we use?
      // Is this an integer type?
      // Can this type generate exponential notation?
      var formatType = formatTypes[type],
          maybeSuffix = /[defgprs%]/.test(type);

      // Set the default precision if not specified,
      // or clamp the specified precision to the supported range.
      // For significant precision, it must be in [1, 21].
      // For fixed precision, it must be in [0, 20].
      precision = precision === undefined ? 6
          : /[gprs]/.test(type) ? Math.max(1, Math.min(21, precision))
          : Math.max(0, Math.min(20, precision));

      function format(value) {
        var valuePrefix = prefix,
            valueSuffix = suffix,
            i, n, c;

        if (type === "c") {
          valueSuffix = formatType(value) + valueSuffix;
          value = "";
        } else {
          value = +value;

          // Determine the sign. -0 is not less than 0, but 1 / -0 is!
          var valueNegative = value < 0 || 1 / value < 0;

          // Perform the initial formatting.
          value = isNaN(value) ? nan : formatType(Math.abs(value), precision);

          // Trim insignificant zeros.
          if (trim) value = formatTrim(value);

          // If a negative value rounds to zero after formatting, and no explicit positive sign is requested, hide the sign.
          if (valueNegative && +value === 0 && sign !== "+") valueNegative = false;

          // Compute the prefix and suffix.
          valuePrefix = (valueNegative ? (sign === "(" ? sign : minus) : sign === "-" || sign === "(" ? "" : sign) + valuePrefix;
          valueSuffix = (type === "s" ? prefixes[8 + prefixExponent / 3] : "") + valueSuffix + (valueNegative && sign === "(" ? ")" : "");

          // Break the formatted value into the integer “value” part that can be
          // grouped, and fractional or exponential “suffix” part that is not.
          if (maybeSuffix) {
            i = -1, n = value.length;
            while (++i < n) {
              if (c = value.charCodeAt(i), 48 > c || c > 57) {
                valueSuffix = (c === 46 ? decimal + value.slice(i + 1) : value.slice(i)) + valueSuffix;
                value = value.slice(0, i);
                break;
              }
            }
          }
        }

        // If the fill character is not "0", grouping is applied before padding.
        if (comma && !zero) value = group(value, Infinity);

        // Compute the padding.
        var length = valuePrefix.length + value.length + valueSuffix.length,
            padding = length < width ? new Array(width - length + 1).join(fill) : "";

        // If the fill character is "0", grouping is applied after padding.
        if (comma && zero) value = group(padding + value, padding.length ? width - valueSuffix.length : Infinity), padding = "";

        // Reconstruct the final output based on the desired alignment.
        switch (align) {
          case "<": value = valuePrefix + value + valueSuffix + padding; break;
          case "=": value = valuePrefix + padding + value + valueSuffix; break;
          case "^": value = padding.slice(0, length = padding.length >> 1) + valuePrefix + value + valueSuffix + padding.slice(length); break;
          default: value = padding + valuePrefix + value + valueSuffix; break;
        }

        return numerals(value);
      }

      format.toString = function() {
        return specifier + "";
      };

      return format;
    }

    function formatPrefix(specifier, value) {
      var f = newFormat((specifier = formatSpecifier(specifier), specifier.type = "f", specifier)),
          e = Math.max(-8, Math.min(8, Math.floor(exponent(value) / 3))) * 3,
          k = Math.pow(10, -e),
          prefix = prefixes[8 + e / 3];
      return function(value) {
        return f(k * value) + prefix;
      };
    }

    return {
      format: newFormat,
      formatPrefix: formatPrefix
    };
  }

  var locale;
  var format;
  var formatPrefix;

  defaultLocale({
    decimal: ".",
    thousands: ",",
    grouping: [3],
    currency: ["$", ""],
    minus: "-"
  });

  function defaultLocale(definition) {
    locale = formatLocale(definition);
    format = locale.format;
    formatPrefix = locale.formatPrefix;
    return locale;
  }

  function precisionFixed(step) {
    return Math.max(0, -exponent(Math.abs(step)));
  }

  function precisionPrefix(step, value) {
    return Math.max(0, Math.max(-8, Math.min(8, Math.floor(exponent(value) / 3))) * 3 - exponent(Math.abs(step)));
  }

  function precisionRound(step, max) {
    step = Math.abs(step), max = Math.abs(max) - step;
    return Math.max(0, exponent(max) - exponent(step)) + 1;
  }

  function tickFormat(start, stop, count, specifier) {
    var step = tickStep(start, stop, count),
        precision;
    specifier = formatSpecifier(specifier == null ? ",f" : specifier);
    switch (specifier.type) {
      case "s": {
        var value = Math.max(Math.abs(start), Math.abs(stop));
        if (specifier.precision == null && !isNaN(precision = precisionPrefix(step, value))) specifier.precision = precision;
        return formatPrefix(specifier, value);
      }
      case "":
      case "e":
      case "g":
      case "p":
      case "r": {
        if (specifier.precision == null && !isNaN(precision = precisionRound(step, Math.max(Math.abs(start), Math.abs(stop))))) specifier.precision = precision - (specifier.type === "e");
        break;
      }
      case "f":
      case "%": {
        if (specifier.precision == null && !isNaN(precision = precisionFixed(step))) specifier.precision = precision - (specifier.type === "%") * 2;
        break;
      }
    }
    return format(specifier);
  }

  function linearish(scale) {
    var domain = scale.domain;

    scale.ticks = function(count) {
      var d = domain();
      return ticks(d[0], d[d.length - 1], count == null ? 10 : count);
    };

    scale.tickFormat = function(count, specifier) {
      var d = domain();
      return tickFormat(d[0], d[d.length - 1], count == null ? 10 : count, specifier);
    };

    scale.nice = function(count) {
      if (count == null) count = 10;

      var d = domain();
      var i0 = 0;
      var i1 = d.length - 1;
      var start = d[i0];
      var stop = d[i1];
      var prestep;
      var step;
      var maxIter = 10;

      if (stop < start) {
        step = start, start = stop, stop = step;
        step = i0, i0 = i1, i1 = step;
      }
      
      while (maxIter-- > 0) {
        step = tickIncrement(start, stop, count);
        if (step === prestep) {
          d[i0] = start;
          d[i1] = stop;
          return domain(d);
        } else if (step > 0) {
          start = Math.floor(start / step) * step;
          stop = Math.ceil(stop / step) * step;
        } else if (step < 0) {
          start = Math.ceil(start * step) / step;
          stop = Math.floor(stop * step) / step;
        } else {
          break;
        }
        prestep = step;
      }

      return scale;
    };

    return scale;
  }

  function linear$1() {
    var scale = continuous();

    scale.copy = function() {
      return copy(scale, linear$1());
    };

    initRange.apply(scale, arguments);

    return linearish(scale);
  }

  var xhtml = "http://www.w3.org/1999/xhtml";

  var namespaces = {
    svg: "http://www.w3.org/2000/svg",
    xhtml: xhtml,
    xlink: "http://www.w3.org/1999/xlink",
    xml: "http://www.w3.org/XML/1998/namespace",
    xmlns: "http://www.w3.org/2000/xmlns/"
  };

  function namespace(name) {
    var prefix = name += "", i = prefix.indexOf(":");
    if (i >= 0 && (prefix = name.slice(0, i)) !== "xmlns") name = name.slice(i + 1);
    return namespaces.hasOwnProperty(prefix) ? {space: namespaces[prefix], local: name} : name; // eslint-disable-line no-prototype-builtins
  }

  function creatorInherit(name) {
    return function() {
      var document = this.ownerDocument,
          uri = this.namespaceURI;
      return uri === xhtml && document.documentElement.namespaceURI === xhtml
          ? document.createElement(name)
          : document.createElementNS(uri, name);
    };
  }

  function creatorFixed(fullname) {
    return function() {
      return this.ownerDocument.createElementNS(fullname.space, fullname.local);
    };
  }

  function creator(name) {
    var fullname = namespace(name);
    return (fullname.local
        ? creatorFixed
        : creatorInherit)(fullname);
  }

  function none() {}

  function selector(selector) {
    return selector == null ? none : function() {
      return this.querySelector(selector);
    };
  }

  function selection_select(select) {
    if (typeof select !== "function") select = selector(select);

    for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, subgroup = subgroups[j] = new Array(n), node, subnode, i = 0; i < n; ++i) {
        if ((node = group[i]) && (subnode = select.call(node, node.__data__, i, group))) {
          if ("__data__" in node) subnode.__data__ = node.__data__;
          subgroup[i] = subnode;
        }
      }
    }

    return new Selection(subgroups, this._parents);
  }

  // Given something array like (or null), returns something that is strictly an
  // array. This is used to ensure that array-like objects passed to d3.selectAll
  // or selection.selectAll are converted into proper arrays when creating a
  // selection; we don’t ever want to create a selection backed by a live
  // HTMLCollection or NodeList. However, note that selection.selectAll will use a
  // static NodeList as a group, since it safely derived from querySelectorAll.
  function array(x) {
    return x == null ? [] : Array.isArray(x) ? x : Array.from(x);
  }

  function empty() {
    return [];
  }

  function selectorAll(selector) {
    return selector == null ? empty : function() {
      return this.querySelectorAll(selector);
    };
  }

  function arrayAll(select) {
    return function() {
      return array(select.apply(this, arguments));
    };
  }

  function selection_selectAll(select) {
    if (typeof select === "function") select = arrayAll(select);
    else select = selectorAll(select);

    for (var groups = this._groups, m = groups.length, subgroups = [], parents = [], j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
        if (node = group[i]) {
          subgroups.push(select.call(node, node.__data__, i, group));
          parents.push(node);
        }
      }
    }

    return new Selection(subgroups, parents);
  }

  function matcher(selector) {
    return function() {
      return this.matches(selector);
    };
  }

  function childMatcher(selector) {
    return function(node) {
      return node.matches(selector);
    };
  }

  var find = Array.prototype.find;

  function childFind(match) {
    return function() {
      return find.call(this.children, match);
    };
  }

  function childFirst() {
    return this.firstElementChild;
  }

  function selection_selectChild(match) {
    return this.select(match == null ? childFirst
        : childFind(typeof match === "function" ? match : childMatcher(match)));
  }

  var filter = Array.prototype.filter;

  function children() {
    return Array.from(this.children);
  }

  function childrenFilter(match) {
    return function() {
      return filter.call(this.children, match);
    };
  }

  function selection_selectChildren(match) {
    return this.selectAll(match == null ? children
        : childrenFilter(typeof match === "function" ? match : childMatcher(match)));
  }

  function selection_filter(match) {
    if (typeof match !== "function") match = matcher(match);

    for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, subgroup = subgroups[j] = [], node, i = 0; i < n; ++i) {
        if ((node = group[i]) && match.call(node, node.__data__, i, group)) {
          subgroup.push(node);
        }
      }
    }

    return new Selection(subgroups, this._parents);
  }

  function sparse(update) {
    return new Array(update.length);
  }

  function selection_enter() {
    return new Selection(this._enter || this._groups.map(sparse), this._parents);
  }

  function EnterNode(parent, datum) {
    this.ownerDocument = parent.ownerDocument;
    this.namespaceURI = parent.namespaceURI;
    this._next = null;
    this._parent = parent;
    this.__data__ = datum;
  }

  EnterNode.prototype = {
    constructor: EnterNode,
    appendChild: function(child) { return this._parent.insertBefore(child, this._next); },
    insertBefore: function(child, next) { return this._parent.insertBefore(child, next); },
    querySelector: function(selector) { return this._parent.querySelector(selector); },
    querySelectorAll: function(selector) { return this._parent.querySelectorAll(selector); }
  };

  function constant$1(x) {
    return function() {
      return x;
    };
  }

  function bindIndex(parent, group, enter, update, exit, data) {
    var i = 0,
        node,
        groupLength = group.length,
        dataLength = data.length;

    // Put any non-null nodes that fit into update.
    // Put any null nodes into enter.
    // Put any remaining data into enter.
    for (; i < dataLength; ++i) {
      if (node = group[i]) {
        node.__data__ = data[i];
        update[i] = node;
      } else {
        enter[i] = new EnterNode(parent, data[i]);
      }
    }

    // Put any non-null nodes that don’t fit into exit.
    for (; i < groupLength; ++i) {
      if (node = group[i]) {
        exit[i] = node;
      }
    }
  }

  function bindKey(parent, group, enter, update, exit, data, key) {
    var i,
        node,
        nodeByKeyValue = new Map,
        groupLength = group.length,
        dataLength = data.length,
        keyValues = new Array(groupLength),
        keyValue;

    // Compute the key for each node.
    // If multiple nodes have the same key, the duplicates are added to exit.
    for (i = 0; i < groupLength; ++i) {
      if (node = group[i]) {
        keyValues[i] = keyValue = key.call(node, node.__data__, i, group) + "";
        if (nodeByKeyValue.has(keyValue)) {
          exit[i] = node;
        } else {
          nodeByKeyValue.set(keyValue, node);
        }
      }
    }

    // Compute the key for each datum.
    // If there a node associated with this key, join and add it to update.
    // If there is not (or the key is a duplicate), add it to enter.
    for (i = 0; i < dataLength; ++i) {
      keyValue = key.call(parent, data[i], i, data) + "";
      if (node = nodeByKeyValue.get(keyValue)) {
        update[i] = node;
        node.__data__ = data[i];
        nodeByKeyValue.delete(keyValue);
      } else {
        enter[i] = new EnterNode(parent, data[i]);
      }
    }

    // Add any remaining nodes that were not bound to data to exit.
    for (i = 0; i < groupLength; ++i) {
      if ((node = group[i]) && (nodeByKeyValue.get(keyValues[i]) === node)) {
        exit[i] = node;
      }
    }
  }

  function datum(node) {
    return node.__data__;
  }

  function selection_data(value, key) {
    if (!arguments.length) return Array.from(this, datum);

    var bind = key ? bindKey : bindIndex,
        parents = this._parents,
        groups = this._groups;

    if (typeof value !== "function") value = constant$1(value);

    for (var m = groups.length, update = new Array(m), enter = new Array(m), exit = new Array(m), j = 0; j < m; ++j) {
      var parent = parents[j],
          group = groups[j],
          groupLength = group.length,
          data = arraylike(value.call(parent, parent && parent.__data__, j, parents)),
          dataLength = data.length,
          enterGroup = enter[j] = new Array(dataLength),
          updateGroup = update[j] = new Array(dataLength),
          exitGroup = exit[j] = new Array(groupLength);

      bind(parent, group, enterGroup, updateGroup, exitGroup, data, key);

      // Now connect the enter nodes to their following update node, such that
      // appendChild can insert the materialized enter node before this node,
      // rather than at the end of the parent node.
      for (var i0 = 0, i1 = 0, previous, next; i0 < dataLength; ++i0) {
        if (previous = enterGroup[i0]) {
          if (i0 >= i1) i1 = i0 + 1;
          while (!(next = updateGroup[i1]) && ++i1 < dataLength);
          previous._next = next || null;
        }
      }
    }

    update = new Selection(update, parents);
    update._enter = enter;
    update._exit = exit;
    return update;
  }

  // Given some data, this returns an array-like view of it: an object that
  // exposes a length property and allows numeric indexing. Note that unlike
  // selectAll, this isn’t worried about “live” collections because the resulting
  // array will only be used briefly while data is being bound. (It is possible to
  // cause the data to change while iterating by using a key function, but please
  // don’t; we’d rather avoid a gratuitous copy.)
  function arraylike(data) {
    return typeof data === "object" && "length" in data
      ? data // Array, TypedArray, NodeList, array-like
      : Array.from(data); // Map, Set, iterable, string, or anything else
  }

  function selection_exit() {
    return new Selection(this._exit || this._groups.map(sparse), this._parents);
  }

  function selection_join(onenter, onupdate, onexit) {
    var enter = this.enter(), update = this, exit = this.exit();
    if (typeof onenter === "function") {
      enter = onenter(enter);
      if (enter) enter = enter.selection();
    } else {
      enter = enter.append(onenter + "");
    }
    if (onupdate != null) {
      update = onupdate(update);
      if (update) update = update.selection();
    }
    if (onexit == null) exit.remove(); else onexit(exit);
    return enter && update ? enter.merge(update).order() : update;
  }

  function selection_merge(context) {
    var selection = context.selection ? context.selection() : context;

    for (var groups0 = this._groups, groups1 = selection._groups, m0 = groups0.length, m1 = groups1.length, m = Math.min(m0, m1), merges = new Array(m0), j = 0; j < m; ++j) {
      for (var group0 = groups0[j], group1 = groups1[j], n = group0.length, merge = merges[j] = new Array(n), node, i = 0; i < n; ++i) {
        if (node = group0[i] || group1[i]) {
          merge[i] = node;
        }
      }
    }

    for (; j < m0; ++j) {
      merges[j] = groups0[j];
    }

    return new Selection(merges, this._parents);
  }

  function selection_order() {

    for (var groups = this._groups, j = -1, m = groups.length; ++j < m;) {
      for (var group = groups[j], i = group.length - 1, next = group[i], node; --i >= 0;) {
        if (node = group[i]) {
          if (next && node.compareDocumentPosition(next) ^ 4) next.parentNode.insertBefore(node, next);
          next = node;
        }
      }
    }

    return this;
  }

  function selection_sort(compare) {
    if (!compare) compare = ascending$1;

    function compareNode(a, b) {
      return a && b ? compare(a.__data__, b.__data__) : !a - !b;
    }

    for (var groups = this._groups, m = groups.length, sortgroups = new Array(m), j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, sortgroup = sortgroups[j] = new Array(n), node, i = 0; i < n; ++i) {
        if (node = group[i]) {
          sortgroup[i] = node;
        }
      }
      sortgroup.sort(compareNode);
    }

    return new Selection(sortgroups, this._parents).order();
  }

  function ascending$1(a, b) {
    return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
  }

  function selection_call() {
    var callback = arguments[0];
    arguments[0] = this;
    callback.apply(null, arguments);
    return this;
  }

  function selection_nodes() {
    return Array.from(this);
  }

  function selection_node() {

    for (var groups = this._groups, j = 0, m = groups.length; j < m; ++j) {
      for (var group = groups[j], i = 0, n = group.length; i < n; ++i) {
        var node = group[i];
        if (node) return node;
      }
    }

    return null;
  }

  function selection_size() {
    let size = 0;
    for (const node of this) ++size; // eslint-disable-line no-unused-vars
    return size;
  }

  function selection_empty() {
    return !this.node();
  }

  function selection_each(callback) {

    for (var groups = this._groups, j = 0, m = groups.length; j < m; ++j) {
      for (var group = groups[j], i = 0, n = group.length, node; i < n; ++i) {
        if (node = group[i]) callback.call(node, node.__data__, i, group);
      }
    }

    return this;
  }

  function attrRemove(name) {
    return function() {
      this.removeAttribute(name);
    };
  }

  function attrRemoveNS(fullname) {
    return function() {
      this.removeAttributeNS(fullname.space, fullname.local);
    };
  }

  function attrConstant(name, value) {
    return function() {
      this.setAttribute(name, value);
    };
  }

  function attrConstantNS(fullname, value) {
    return function() {
      this.setAttributeNS(fullname.space, fullname.local, value);
    };
  }

  function attrFunction(name, value) {
    return function() {
      var v = value.apply(this, arguments);
      if (v == null) this.removeAttribute(name);
      else this.setAttribute(name, v);
    };
  }

  function attrFunctionNS(fullname, value) {
    return function() {
      var v = value.apply(this, arguments);
      if (v == null) this.removeAttributeNS(fullname.space, fullname.local);
      else this.setAttributeNS(fullname.space, fullname.local, v);
    };
  }

  function selection_attr(name, value) {
    var fullname = namespace(name);

    if (arguments.length < 2) {
      var node = this.node();
      return fullname.local
          ? node.getAttributeNS(fullname.space, fullname.local)
          : node.getAttribute(fullname);
    }

    return this.each((value == null
        ? (fullname.local ? attrRemoveNS : attrRemove) : (typeof value === "function"
        ? (fullname.local ? attrFunctionNS : attrFunction)
        : (fullname.local ? attrConstantNS : attrConstant)))(fullname, value));
  }

  function defaultView(node) {
    return (node.ownerDocument && node.ownerDocument.defaultView) // node is a Node
        || (node.document && node) // node is a Window
        || node.defaultView; // node is a Document
  }

  function styleRemove(name) {
    return function() {
      this.style.removeProperty(name);
    };
  }

  function styleConstant(name, value, priority) {
    return function() {
      this.style.setProperty(name, value, priority);
    };
  }

  function styleFunction(name, value, priority) {
    return function() {
      var v = value.apply(this, arguments);
      if (v == null) this.style.removeProperty(name);
      else this.style.setProperty(name, v, priority);
    };
  }

  function selection_style(name, value, priority) {
    return arguments.length > 1
        ? this.each((value == null
              ? styleRemove : typeof value === "function"
              ? styleFunction
              : styleConstant)(name, value, priority == null ? "" : priority))
        : styleValue(this.node(), name);
  }

  function styleValue(node, name) {
    return node.style.getPropertyValue(name)
        || defaultView(node).getComputedStyle(node, null).getPropertyValue(name);
  }

  function propertyRemove(name) {
    return function() {
      delete this[name];
    };
  }

  function propertyConstant(name, value) {
    return function() {
      this[name] = value;
    };
  }

  function propertyFunction(name, value) {
    return function() {
      var v = value.apply(this, arguments);
      if (v == null) delete this[name];
      else this[name] = v;
    };
  }

  function selection_property(name, value) {
    return arguments.length > 1
        ? this.each((value == null
            ? propertyRemove : typeof value === "function"
            ? propertyFunction
            : propertyConstant)(name, value))
        : this.node()[name];
  }

  function classArray(string) {
    return string.trim().split(/^|\s+/);
  }

  function classList(node) {
    return node.classList || new ClassList(node);
  }

  function ClassList(node) {
    this._node = node;
    this._names = classArray(node.getAttribute("class") || "");
  }

  ClassList.prototype = {
    add: function(name) {
      var i = this._names.indexOf(name);
      if (i < 0) {
        this._names.push(name);
        this._node.setAttribute("class", this._names.join(" "));
      }
    },
    remove: function(name) {
      var i = this._names.indexOf(name);
      if (i >= 0) {
        this._names.splice(i, 1);
        this._node.setAttribute("class", this._names.join(" "));
      }
    },
    contains: function(name) {
      return this._names.indexOf(name) >= 0;
    }
  };

  function classedAdd(node, names) {
    var list = classList(node), i = -1, n = names.length;
    while (++i < n) list.add(names[i]);
  }

  function classedRemove(node, names) {
    var list = classList(node), i = -1, n = names.length;
    while (++i < n) list.remove(names[i]);
  }

  function classedTrue(names) {
    return function() {
      classedAdd(this, names);
    };
  }

  function classedFalse(names) {
    return function() {
      classedRemove(this, names);
    };
  }

  function classedFunction(names, value) {
    return function() {
      (value.apply(this, arguments) ? classedAdd : classedRemove)(this, names);
    };
  }

  function selection_classed(name, value) {
    var names = classArray(name + "");

    if (arguments.length < 2) {
      var list = classList(this.node()), i = -1, n = names.length;
      while (++i < n) if (!list.contains(names[i])) return false;
      return true;
    }

    return this.each((typeof value === "function"
        ? classedFunction : value
        ? classedTrue
        : classedFalse)(names, value));
  }

  function textRemove() {
    this.textContent = "";
  }

  function textConstant(value) {
    return function() {
      this.textContent = value;
    };
  }

  function textFunction(value) {
    return function() {
      var v = value.apply(this, arguments);
      this.textContent = v == null ? "" : v;
    };
  }

  function selection_text(value) {
    return arguments.length
        ? this.each(value == null
            ? textRemove : (typeof value === "function"
            ? textFunction
            : textConstant)(value))
        : this.node().textContent;
  }

  function htmlRemove() {
    this.innerHTML = "";
  }

  function htmlConstant(value) {
    return function() {
      this.innerHTML = value;
    };
  }

  function htmlFunction(value) {
    return function() {
      var v = value.apply(this, arguments);
      this.innerHTML = v == null ? "" : v;
    };
  }

  function selection_html(value) {
    return arguments.length
        ? this.each(value == null
            ? htmlRemove : (typeof value === "function"
            ? htmlFunction
            : htmlConstant)(value))
        : this.node().innerHTML;
  }

  function raise() {
    if (this.nextSibling) this.parentNode.appendChild(this);
  }

  function selection_raise() {
    return this.each(raise);
  }

  function lower() {
    if (this.previousSibling) this.parentNode.insertBefore(this, this.parentNode.firstChild);
  }

  function selection_lower() {
    return this.each(lower);
  }

  function selection_append(name) {
    var create = typeof name === "function" ? name : creator(name);
    return this.select(function() {
      return this.appendChild(create.apply(this, arguments));
    });
  }

  function constantNull() {
    return null;
  }

  function selection_insert(name, before) {
    var create = typeof name === "function" ? name : creator(name),
        select = before == null ? constantNull : typeof before === "function" ? before : selector(before);
    return this.select(function() {
      return this.insertBefore(create.apply(this, arguments), select.apply(this, arguments) || null);
    });
  }

  function remove() {
    var parent = this.parentNode;
    if (parent) parent.removeChild(this);
  }

  function selection_remove() {
    return this.each(remove);
  }

  function selection_cloneShallow() {
    var clone = this.cloneNode(false), parent = this.parentNode;
    return parent ? parent.insertBefore(clone, this.nextSibling) : clone;
  }

  function selection_cloneDeep() {
    var clone = this.cloneNode(true), parent = this.parentNode;
    return parent ? parent.insertBefore(clone, this.nextSibling) : clone;
  }

  function selection_clone(deep) {
    return this.select(deep ? selection_cloneDeep : selection_cloneShallow);
  }

  function selection_datum(value) {
    return arguments.length
        ? this.property("__data__", value)
        : this.node().__data__;
  }

  function contextListener(listener) {
    return function(event) {
      listener.call(this, event, this.__data__);
    };
  }

  function parseTypenames(typenames) {
    return typenames.trim().split(/^|\s+/).map(function(t) {
      var name = "", i = t.indexOf(".");
      if (i >= 0) name = t.slice(i + 1), t = t.slice(0, i);
      return {type: t, name: name};
    });
  }

  function onRemove(typename) {
    return function() {
      var on = this.__on;
      if (!on) return;
      for (var j = 0, i = -1, m = on.length, o; j < m; ++j) {
        if (o = on[j], (!typename.type || o.type === typename.type) && o.name === typename.name) {
          this.removeEventListener(o.type, o.listener, o.options);
        } else {
          on[++i] = o;
        }
      }
      if (++i) on.length = i;
      else delete this.__on;
    };
  }

  function onAdd(typename, value, options) {
    return function() {
      var on = this.__on, o, listener = contextListener(value);
      if (on) for (var j = 0, m = on.length; j < m; ++j) {
        if ((o = on[j]).type === typename.type && o.name === typename.name) {
          this.removeEventListener(o.type, o.listener, o.options);
          this.addEventListener(o.type, o.listener = listener, o.options = options);
          o.value = value;
          return;
        }
      }
      this.addEventListener(typename.type, listener, options);
      o = {type: typename.type, name: typename.name, value: value, listener: listener, options: options};
      if (!on) this.__on = [o];
      else on.push(o);
    };
  }

  function selection_on(typename, value, options) {
    var typenames = parseTypenames(typename + ""), i, n = typenames.length, t;

    if (arguments.length < 2) {
      var on = this.node().__on;
      if (on) for (var j = 0, m = on.length, o; j < m; ++j) {
        for (i = 0, o = on[j]; i < n; ++i) {
          if ((t = typenames[i]).type === o.type && t.name === o.name) {
            return o.value;
          }
        }
      }
      return;
    }

    on = value ? onAdd : onRemove;
    for (i = 0; i < n; ++i) this.each(on(typenames[i], value, options));
    return this;
  }

  function dispatchEvent(node, type, params) {
    var window = defaultView(node),
        event = window.CustomEvent;

    if (typeof event === "function") {
      event = new event(type, params);
    } else {
      event = window.document.createEvent("Event");
      if (params) event.initEvent(type, params.bubbles, params.cancelable), event.detail = params.detail;
      else event.initEvent(type, false, false);
    }

    node.dispatchEvent(event);
  }

  function dispatchConstant(type, params) {
    return function() {
      return dispatchEvent(this, type, params);
    };
  }

  function dispatchFunction(type, params) {
    return function() {
      return dispatchEvent(this, type, params.apply(this, arguments));
    };
  }

  function selection_dispatch(type, params) {
    return this.each((typeof params === "function"
        ? dispatchFunction
        : dispatchConstant)(type, params));
  }

  function* selection_iterator() {
    for (var groups = this._groups, j = 0, m = groups.length; j < m; ++j) {
      for (var group = groups[j], i = 0, n = group.length, node; i < n; ++i) {
        if (node = group[i]) yield node;
      }
    }
  }

  var root = [null];

  function Selection(groups, parents) {
    this._groups = groups;
    this._parents = parents;
  }

  function selection() {
    return new Selection([[document.documentElement]], root);
  }

  function selection_selection() {
    return this;
  }

  Selection.prototype = selection.prototype = {
    constructor: Selection,
    select: selection_select,
    selectAll: selection_selectAll,
    selectChild: selection_selectChild,
    selectChildren: selection_selectChildren,
    filter: selection_filter,
    data: selection_data,
    enter: selection_enter,
    exit: selection_exit,
    join: selection_join,
    merge: selection_merge,
    selection: selection_selection,
    order: selection_order,
    sort: selection_sort,
    call: selection_call,
    nodes: selection_nodes,
    node: selection_node,
    size: selection_size,
    empty: selection_empty,
    each: selection_each,
    attr: selection_attr,
    style: selection_style,
    property: selection_property,
    classed: selection_classed,
    text: selection_text,
    html: selection_html,
    raise: selection_raise,
    lower: selection_lower,
    append: selection_append,
    insert: selection_insert,
    remove: selection_remove,
    clone: selection_clone,
    datum: selection_datum,
    on: selection_on,
    dispatch: selection_dispatch,
    [Symbol.iterator]: selection_iterator
  };

  function select(selector) {
    return typeof selector === "string"
        ? new Selection([[document.querySelector(selector)]], [document.documentElement])
        : new Selection([[selector]], root);
  }

  function sourceEvent(event) {
    let sourceEvent;
    while (sourceEvent = event.sourceEvent) event = sourceEvent;
    return event;
  }

  function pointer(event, node) {
    event = sourceEvent(event);
    if (node === undefined) node = event.currentTarget;
    if (node) {
      var svg = node.ownerSVGElement || node;
      if (svg.createSVGPoint) {
        var point = svg.createSVGPoint();
        point.x = event.clientX, point.y = event.clientY;
        point = point.matrixTransform(node.getScreenCTM().inverse());
        return [point.x, point.y];
      }
      if (node.getBoundingClientRect) {
        var rect = node.getBoundingClientRect();
        return [event.clientX - rect.left - node.clientLeft, event.clientY - rect.top - node.clientTop];
      }
    }
    return [event.pageX, event.pageY];
  }

  const pi = Math.PI,
      tau = 2 * pi,
      epsilon$1 = 1e-6,
      tauEpsilon = tau - epsilon$1;

  function Path() {
    this._x0 = this._y0 = // start of current subpath
    this._x1 = this._y1 = null; // end of current subpath
    this._ = "";
  }

  function path() {
    return new Path;
  }

  Path.prototype = path.prototype = {
    constructor: Path,
    moveTo: function(x, y) {
      this._ += "M" + (this._x0 = this._x1 = +x) + "," + (this._y0 = this._y1 = +y);
    },
    closePath: function() {
      if (this._x1 !== null) {
        this._x1 = this._x0, this._y1 = this._y0;
        this._ += "Z";
      }
    },
    lineTo: function(x, y) {
      this._ += "L" + (this._x1 = +x) + "," + (this._y1 = +y);
    },
    quadraticCurveTo: function(x1, y1, x, y) {
      this._ += "Q" + (+x1) + "," + (+y1) + "," + (this._x1 = +x) + "," + (this._y1 = +y);
    },
    bezierCurveTo: function(x1, y1, x2, y2, x, y) {
      this._ += "C" + (+x1) + "," + (+y1) + "," + (+x2) + "," + (+y2) + "," + (this._x1 = +x) + "," + (this._y1 = +y);
    },
    arcTo: function(x1, y1, x2, y2, r) {
      x1 = +x1, y1 = +y1, x2 = +x2, y2 = +y2, r = +r;
      var x0 = this._x1,
          y0 = this._y1,
          x21 = x2 - x1,
          y21 = y2 - y1,
          x01 = x0 - x1,
          y01 = y0 - y1,
          l01_2 = x01 * x01 + y01 * y01;

      // Is the radius negative? Error.
      if (r < 0) throw new Error("negative radius: " + r);

      // Is this path empty? Move to (x1,y1).
      if (this._x1 === null) {
        this._ += "M" + (this._x1 = x1) + "," + (this._y1 = y1);
      }

      // Or, is (x1,y1) coincident with (x0,y0)? Do nothing.
      else if (!(l01_2 > epsilon$1));

      // Or, are (x0,y0), (x1,y1) and (x2,y2) collinear?
      // Equivalently, is (x1,y1) coincident with (x2,y2)?
      // Or, is the radius zero? Line to (x1,y1).
      else if (!(Math.abs(y01 * x21 - y21 * x01) > epsilon$1) || !r) {
        this._ += "L" + (this._x1 = x1) + "," + (this._y1 = y1);
      }

      // Otherwise, draw an arc!
      else {
        var x20 = x2 - x0,
            y20 = y2 - y0,
            l21_2 = x21 * x21 + y21 * y21,
            l20_2 = x20 * x20 + y20 * y20,
            l21 = Math.sqrt(l21_2),
            l01 = Math.sqrt(l01_2),
            l = r * Math.tan((pi - Math.acos((l21_2 + l01_2 - l20_2) / (2 * l21 * l01))) / 2),
            t01 = l / l01,
            t21 = l / l21;

        // If the start tangent is not coincident with (x0,y0), line to.
        if (Math.abs(t01 - 1) > epsilon$1) {
          this._ += "L" + (x1 + t01 * x01) + "," + (y1 + t01 * y01);
        }

        this._ += "A" + r + "," + r + ",0,0," + (+(y01 * x20 > x01 * y20)) + "," + (this._x1 = x1 + t21 * x21) + "," + (this._y1 = y1 + t21 * y21);
      }
    },
    arc: function(x, y, r, a0, a1, ccw) {
      x = +x, y = +y, r = +r, ccw = !!ccw;
      var dx = r * Math.cos(a0),
          dy = r * Math.sin(a0),
          x0 = x + dx,
          y0 = y + dy,
          cw = 1 ^ ccw,
          da = ccw ? a0 - a1 : a1 - a0;

      // Is the radius negative? Error.
      if (r < 0) throw new Error("negative radius: " + r);

      // Is this path empty? Move to (x0,y0).
      if (this._x1 === null) {
        this._ += "M" + x0 + "," + y0;
      }

      // Or, is (x0,y0) not coincident with the previous point? Line to (x0,y0).
      else if (Math.abs(this._x1 - x0) > epsilon$1 || Math.abs(this._y1 - y0) > epsilon$1) {
        this._ += "L" + x0 + "," + y0;
      }

      // Is this arc empty? We’re done.
      if (!r) return;

      // Does the angle go the wrong way? Flip the direction.
      if (da < 0) da = da % tau + tau;

      // Is this a complete circle? Draw two arcs to complete the circle.
      if (da > tauEpsilon) {
        this._ += "A" + r + "," + r + ",0,1," + cw + "," + (x - dx) + "," + (y - dy) + "A" + r + "," + r + ",0,1," + cw + "," + (this._x1 = x0) + "," + (this._y1 = y0);
      }

      // Is this arc non-empty? Draw an arc!
      else if (da > epsilon$1) {
        this._ += "A" + r + "," + r + ",0," + (+(da >= pi)) + "," + cw + "," + (this._x1 = x + r * Math.cos(a1)) + "," + (this._y1 = y + r * Math.sin(a1));
      }
    },
    rect: function(x, y, w, h) {
      this._ += "M" + (this._x0 = this._x1 = +x) + "," + (this._y0 = this._y1 = +y) + "h" + (+w) + "v" + (+h) + "h" + (-w) + "Z";
    },
    toString: function() {
      return this._;
    }
  };

  function constant$2(x) {
    return function constant() {
      return x;
    };
  }

  var slice = Array.prototype.slice;

  function x(p) {
    return p[0];
  }

  function y(p) {
    return p[1];
  }

  function pointRadial(x, y) {
    return [(y = +y) * Math.cos(x -= Math.PI / 2), y * Math.sin(x)];
  }

  class Bump {
    constructor(context, x) {
      this._context = context;
      this._x = x;
    }
    areaStart() {
      this._line = 0;
    }
    areaEnd() {
      this._line = NaN;
    }
    lineStart() {
      this._point = 0;
    }
    lineEnd() {
      if (this._line || (this._line !== 0 && this._point === 1)) this._context.closePath();
      this._line = 1 - this._line;
    }
    point(x, y) {
      x = +x, y = +y;
      switch (this._point) {
        case 0: {
          this._point = 1;
          if (this._line) this._context.lineTo(x, y);
          else this._context.moveTo(x, y);
          break;
        }
        case 1: this._point = 2; // falls through
        default: {
          if (this._x) this._context.bezierCurveTo(this._x0 = (this._x0 + x) / 2, this._y0, this._x0, y, x, y);
          else this._context.bezierCurveTo(this._x0, this._y0 = (this._y0 + y) / 2, x, this._y0, x, y);
          break;
        }
      }
      this._x0 = x, this._y0 = y;
    }
  }

  class BumpRadial {
    constructor(context) {
      this._context = context;
    }
    lineStart() {
      this._point = 0;
    }
    lineEnd() {}
    point(x, y) {
      x = +x, y = +y;
      if (this._point++ === 0) {
        this._x0 = x, this._y0 = y;
      } else {
        const p0 = pointRadial(this._x0, this._y0);
        const p1 = pointRadial(this._x0, this._y0 = (this._y0 + y) / 2);
        const p2 = pointRadial(x, this._y0);
        const p3 = pointRadial(x, y);
        this._context.moveTo(...p0);
        this._context.bezierCurveTo(...p1, ...p2, ...p3);
      }
    }
  }

  function bumpX(context) {
    return new Bump(context, true);
  }

  function bumpY(context) {
    return new Bump(context, false);
  }

  function bumpRadial(context) {
    return new BumpRadial(context);
  }

  function linkSource(d) {
    return d.source;
  }

  function linkTarget(d) {
    return d.target;
  }

  function link(curve) {
    let source = linkSource;
    let target = linkTarget;
    let x$1 = x;
    let y$1 = y;
    let context = null;
    let output = null;

    function link() {
      let buffer;
      const argv = slice.call(arguments);
      const s = source.apply(this, argv);
      const t = target.apply(this, argv);
      if (context == null) output = curve(buffer = path());
      output.lineStart();
      argv[0] = s, output.point(+x$1.apply(this, argv), +y$1.apply(this, argv));
      argv[0] = t, output.point(+x$1.apply(this, argv), +y$1.apply(this, argv));
      output.lineEnd();
      if (buffer) return output = null, buffer + "" || null;
    }

    link.source = function(_) {
      return arguments.length ? (source = _, link) : source;
    };

    link.target = function(_) {
      return arguments.length ? (target = _, link) : target;
    };

    link.x = function(_) {
      return arguments.length ? (x$1 = typeof _ === "function" ? _ : constant$2(+_), link) : x$1;
    };

    link.y = function(_) {
      return arguments.length ? (y$1 = typeof _ === "function" ? _ : constant$2(+_), link) : y$1;
    };

    link.context = function(_) {
      return arguments.length ? (_ == null ? context = output = null : output = curve(context = _), link) : context;
    };

    return link;
  }

  function linkHorizontal() {
    return link(bumpX);
  }

  function linkVertical() {
    return link(bumpY);
  }

  function linkRadial() {
    const l = link(bumpRadial);
    l.angle = l.x, delete l.x;
    l.radius = l.y, delete l.y;
    return l;
  }

  var noop = {value: function() {}};

  function dispatch() {
    for (var i = 0, n = arguments.length, _ = {}, t; i < n; ++i) {
      if (!(t = arguments[i] + "") || (t in _) || /[\s.]/.test(t)) throw new Error("illegal type: " + t);
      _[t] = [];
    }
    return new Dispatch(_);
  }

  function Dispatch(_) {
    this._ = _;
  }

  function parseTypenames$1(typenames, types) {
    return typenames.trim().split(/^|\s+/).map(function(t) {
      var name = "", i = t.indexOf(".");
      if (i >= 0) name = t.slice(i + 1), t = t.slice(0, i);
      if (t && !types.hasOwnProperty(t)) throw new Error("unknown type: " + t);
      return {type: t, name: name};
    });
  }

  Dispatch.prototype = dispatch.prototype = {
    constructor: Dispatch,
    on: function(typename, callback) {
      var _ = this._,
          T = parseTypenames$1(typename + "", _),
          t,
          i = -1,
          n = T.length;

      // If no callback was specified, return the callback of the given type and name.
      if (arguments.length < 2) {
        while (++i < n) if ((t = (typename = T[i]).type) && (t = get(_[t], typename.name))) return t;
        return;
      }

      // If a type was specified, set the callback for the given type and name.
      // Otherwise, if a null callback was specified, remove callbacks of the given name.
      if (callback != null && typeof callback !== "function") throw new Error("invalid callback: " + callback);
      while (++i < n) {
        if (t = (typename = T[i]).type) _[t] = set(_[t], typename.name, callback);
        else if (callback == null) for (t in _) _[t] = set(_[t], typename.name, null);
      }

      return this;
    },
    copy: function() {
      var copy = {}, _ = this._;
      for (var t in _) copy[t] = _[t].slice();
      return new Dispatch(copy);
    },
    call: function(type, that) {
      if ((n = arguments.length - 2) > 0) for (var args = new Array(n), i = 0, n, t; i < n; ++i) args[i] = arguments[i + 2];
      if (!this._.hasOwnProperty(type)) throw new Error("unknown type: " + type);
      for (t = this._[type], i = 0, n = t.length; i < n; ++i) t[i].value.apply(that, args);
    },
    apply: function(type, that, args) {
      if (!this._.hasOwnProperty(type)) throw new Error("unknown type: " + type);
      for (var t = this._[type], i = 0, n = t.length; i < n; ++i) t[i].value.apply(that, args);
    }
  };

  function get(type, name) {
    for (var i = 0, n = type.length, c; i < n; ++i) {
      if ((c = type[i]).name === name) {
        return c.value;
      }
    }
  }

  function set(type, name, callback) {
    for (var i = 0, n = type.length; i < n; ++i) {
      if (type[i].name === name) {
        type[i] = noop, type = type.slice(0, i).concat(type.slice(i + 1));
        break;
      }
    }
    if (callback != null) type.push({name: name, value: callback});
    return type;
  }

  // These are typically used in conjunction with noevent to ensure that we can
  const nonpassivecapture = {capture: true, passive: false};

  function noevent(event) {
    event.preventDefault();
    event.stopImmediatePropagation();
  }

  function dragDisable(view) {
    var root = view.document.documentElement,
        selection = select(view).on("dragstart.drag", noevent, nonpassivecapture);
    if ("onselectstart" in root) {
      selection.on("selectstart.drag", noevent, nonpassivecapture);
    } else {
      root.__noselect = root.style.MozUserSelect;
      root.style.MozUserSelect = "none";
    }
  }

  function yesdrag(view, noclick) {
    var root = view.document.documentElement,
        selection = select(view).on("dragstart.drag", null);
    if (noclick) {
      selection.on("click.drag", noevent, nonpassivecapture);
      setTimeout(function() { selection.on("click.drag", null); }, 0);
    }
    if ("onselectstart" in root) {
      selection.on("selectstart.drag", null);
    } else {
      root.style.MozUserSelect = root.__noselect;
      delete root.__noselect;
    }
  }

  var frame = 0, // is an animation frame pending?
      timeout = 0, // is a timeout pending?
      interval = 0, // are any timers active?
      pokeDelay = 1000, // how frequently we check for clock skew
      taskHead,
      taskTail,
      clockLast = 0,
      clockNow = 0,
      clockSkew = 0,
      clock = typeof performance === "object" && performance.now ? performance : Date,
      setFrame = typeof window === "object" && window.requestAnimationFrame ? window.requestAnimationFrame.bind(window) : function(f) { setTimeout(f, 17); };

  function now() {
    return clockNow || (setFrame(clearNow), clockNow = clock.now() + clockSkew);
  }

  function clearNow() {
    clockNow = 0;
  }

  function Timer() {
    this._call =
    this._time =
    this._next = null;
  }

  Timer.prototype = timer.prototype = {
    constructor: Timer,
    restart: function(callback, delay, time) {
      if (typeof callback !== "function") throw new TypeError("callback is not a function");
      time = (time == null ? now() : +time) + (delay == null ? 0 : +delay);
      if (!this._next && taskTail !== this) {
        if (taskTail) taskTail._next = this;
        else taskHead = this;
        taskTail = this;
      }
      this._call = callback;
      this._time = time;
      sleep();
    },
    stop: function() {
      if (this._call) {
        this._call = null;
        this._time = Infinity;
        sleep();
      }
    }
  };

  function timer(callback, delay, time) {
    var t = new Timer;
    t.restart(callback, delay, time);
    return t;
  }

  function timerFlush() {
    now(); // Get the current time, if not already set.
    ++frame; // Pretend we’ve set an alarm, if we haven’t already.
    var t = taskHead, e;
    while (t) {
      if ((e = clockNow - t._time) >= 0) t._call.call(null, e);
      t = t._next;
    }
    --frame;
  }

  function wake() {
    clockNow = (clockLast = clock.now()) + clockSkew;
    frame = timeout = 0;
    try {
      timerFlush();
    } finally {
      frame = 0;
      nap();
      clockNow = 0;
    }
  }

  function poke() {
    var now = clock.now(), delay = now - clockLast;
    if (delay > pokeDelay) clockSkew -= delay, clockLast = now;
  }

  function nap() {
    var t0, t1 = taskHead, t2, time = Infinity;
    while (t1) {
      if (t1._call) {
        if (time > t1._time) time = t1._time;
        t0 = t1, t1 = t1._next;
      } else {
        t2 = t1._next, t1._next = null;
        t1 = t0 ? t0._next = t2 : taskHead = t2;
      }
    }
    taskTail = t0;
    sleep(time);
  }

  function sleep(time) {
    if (frame) return; // Soonest alarm already set, or will be.
    if (timeout) timeout = clearTimeout(timeout);
    var delay = time - clockNow; // Strictly less than if we recomputed clockNow.
    if (delay > 24) {
      if (time < Infinity) timeout = setTimeout(wake, time - clock.now() - clockSkew);
      if (interval) interval = clearInterval(interval);
    } else {
      if (!interval) clockLast = clock.now(), interval = setInterval(poke, pokeDelay);
      frame = 1, setFrame(wake);
    }
  }

  function timeout$1(callback, delay, time) {
    var t = new Timer;
    delay = delay == null ? 0 : +delay;
    t.restart(function(elapsed) {
      t.stop();
      callback(elapsed + delay);
    }, delay, time);
    return t;
  }

  var emptyOn = dispatch("start", "end", "cancel", "interrupt");
  var emptyTween = [];

  var CREATED = 0;
  var SCHEDULED = 1;
  var STARTING = 2;
  var STARTED = 3;
  var RUNNING = 4;
  var ENDING = 5;
  var ENDED = 6;

  function schedule(node, name, id, index, group, timing) {
    var schedules = node.__transition;
    if (!schedules) node.__transition = {};
    else if (id in schedules) return;
    create(node, id, {
      name: name,
      index: index, // For context during callback.
      group: group, // For context during callback.
      on: emptyOn,
      tween: emptyTween,
      time: timing.time,
      delay: timing.delay,
      duration: timing.duration,
      ease: timing.ease,
      timer: null,
      state: CREATED
    });
  }

  function init(node, id) {
    var schedule = get$1(node, id);
    if (schedule.state > CREATED) throw new Error("too late; already scheduled");
    return schedule;
  }

  function set$1(node, id) {
    var schedule = get$1(node, id);
    if (schedule.state > STARTED) throw new Error("too late; already running");
    return schedule;
  }

  function get$1(node, id) {
    var schedule = node.__transition;
    if (!schedule || !(schedule = schedule[id])) throw new Error("transition not found");
    return schedule;
  }

  function create(node, id, self) {
    var schedules = node.__transition,
        tween;

    // Initialize the self timer when the transition is created.
    // Note the actual delay is not known until the first callback!
    schedules[id] = self;
    self.timer = timer(schedule, 0, self.time);

    function schedule(elapsed) {
      self.state = SCHEDULED;
      self.timer.restart(start, self.delay, self.time);

      // If the elapsed delay is less than our first sleep, start immediately.
      if (self.delay <= elapsed) start(elapsed - self.delay);
    }

    function start(elapsed) {
      var i, j, n, o;

      // If the state is not SCHEDULED, then we previously errored on start.
      if (self.state !== SCHEDULED) return stop();

      for (i in schedules) {
        o = schedules[i];
        if (o.name !== self.name) continue;

        // While this element already has a starting transition during this frame,
        // defer starting an interrupting transition until that transition has a
        // chance to tick (and possibly end); see d3/d3-transition#54!
        if (o.state === STARTED) return timeout$1(start);

        // Interrupt the active transition, if any.
        if (o.state === RUNNING) {
          o.state = ENDED;
          o.timer.stop();
          o.on.call("interrupt", node, node.__data__, o.index, o.group);
          delete schedules[i];
        }

        // Cancel any pre-empted transitions.
        else if (+i < id) {
          o.state = ENDED;
          o.timer.stop();
          o.on.call("cancel", node, node.__data__, o.index, o.group);
          delete schedules[i];
        }
      }

      // Defer the first tick to end of the current frame; see d3/d3#1576.
      // Note the transition may be canceled after start and before the first tick!
      // Note this must be scheduled before the start event; see d3/d3-transition#16!
      // Assuming this is successful, subsequent callbacks go straight to tick.
      timeout$1(function() {
        if (self.state === STARTED) {
          self.state = RUNNING;
          self.timer.restart(tick, self.delay, self.time);
          tick(elapsed);
        }
      });

      // Dispatch the start event.
      // Note this must be done before the tween are initialized.
      self.state = STARTING;
      self.on.call("start", node, node.__data__, self.index, self.group);
      if (self.state !== STARTING) return; // interrupted
      self.state = STARTED;

      // Initialize the tween, deleting null tween.
      tween = new Array(n = self.tween.length);
      for (i = 0, j = -1; i < n; ++i) {
        if (o = self.tween[i].value.call(node, node.__data__, self.index, self.group)) {
          tween[++j] = o;
        }
      }
      tween.length = j + 1;
    }

    function tick(elapsed) {
      var t = elapsed < self.duration ? self.ease.call(null, elapsed / self.duration) : (self.timer.restart(stop), self.state = ENDING, 1),
          i = -1,
          n = tween.length;

      while (++i < n) {
        tween[i].call(node, t);
      }

      // Dispatch the end event.
      if (self.state === ENDING) {
        self.on.call("end", node, node.__data__, self.index, self.group);
        stop();
      }
    }

    function stop() {
      self.state = ENDED;
      self.timer.stop();
      delete schedules[id];
      for (var i in schedules) return; // eslint-disable-line no-unused-vars
      delete node.__transition;
    }
  }

  function interrupt(node, name) {
    var schedules = node.__transition,
        schedule,
        active,
        empty = true,
        i;

    if (!schedules) return;

    name = name == null ? null : name + "";

    for (i in schedules) {
      if ((schedule = schedules[i]).name !== name) { empty = false; continue; }
      active = schedule.state > STARTING && schedule.state < ENDING;
      schedule.state = ENDED;
      schedule.timer.stop();
      schedule.on.call(active ? "interrupt" : "cancel", node, node.__data__, schedule.index, schedule.group);
      delete schedules[i];
    }

    if (empty) delete node.__transition;
  }

  function selection_interrupt(name) {
    return this.each(function() {
      interrupt(this, name);
    });
  }

  function tweenRemove(id, name) {
    var tween0, tween1;
    return function() {
      var schedule = set$1(this, id),
          tween = schedule.tween;

      // If this node shared tween with the previous node,
      // just assign the updated shared tween and we’re done!
      // Otherwise, copy-on-write.
      if (tween !== tween0) {
        tween1 = tween0 = tween;
        for (var i = 0, n = tween1.length; i < n; ++i) {
          if (tween1[i].name === name) {
            tween1 = tween1.slice();
            tween1.splice(i, 1);
            break;
          }
        }
      }

      schedule.tween = tween1;
    };
  }

  function tweenFunction(id, name, value) {
    var tween0, tween1;
    if (typeof value !== "function") throw new Error;
    return function() {
      var schedule = set$1(this, id),
          tween = schedule.tween;

      // If this node shared tween with the previous node,
      // just assign the updated shared tween and we’re done!
      // Otherwise, copy-on-write.
      if (tween !== tween0) {
        tween1 = (tween0 = tween).slice();
        for (var t = {name: name, value: value}, i = 0, n = tween1.length; i < n; ++i) {
          if (tween1[i].name === name) {
            tween1[i] = t;
            break;
          }
        }
        if (i === n) tween1.push(t);
      }

      schedule.tween = tween1;
    };
  }

  function transition_tween(name, value) {
    var id = this._id;

    name += "";

    if (arguments.length < 2) {
      var tween = get$1(this.node(), id).tween;
      for (var i = 0, n = tween.length, t; i < n; ++i) {
        if ((t = tween[i]).name === name) {
          return t.value;
        }
      }
      return null;
    }

    return this.each((value == null ? tweenRemove : tweenFunction)(id, name, value));
  }

  function tweenValue(transition, name, value) {
    var id = transition._id;

    transition.each(function() {
      var schedule = set$1(this, id);
      (schedule.value || (schedule.value = {}))[name] = value.apply(this, arguments);
    });

    return function(node) {
      return get$1(node, id).value[name];
    };
  }

  function interpolate$1(a, b) {
    var c;
    return (typeof b === "number" ? interpolateNumber
        : b instanceof color ? interpolateRgb
        : (c = color(b)) ? (b = c, interpolateRgb)
        : interpolateString)(a, b);
  }

  function attrRemove$1(name) {
    return function() {
      this.removeAttribute(name);
    };
  }

  function attrRemoveNS$1(fullname) {
    return function() {
      this.removeAttributeNS(fullname.space, fullname.local);
    };
  }

  function attrConstant$1(name, interpolate, value1) {
    var string00,
        string1 = value1 + "",
        interpolate0;
    return function() {
      var string0 = this.getAttribute(name);
      return string0 === string1 ? null
          : string0 === string00 ? interpolate0
          : interpolate0 = interpolate(string00 = string0, value1);
    };
  }

  function attrConstantNS$1(fullname, interpolate, value1) {
    var string00,
        string1 = value1 + "",
        interpolate0;
    return function() {
      var string0 = this.getAttributeNS(fullname.space, fullname.local);
      return string0 === string1 ? null
          : string0 === string00 ? interpolate0
          : interpolate0 = interpolate(string00 = string0, value1);
    };
  }

  function attrFunction$1(name, interpolate, value) {
    var string00,
        string10,
        interpolate0;
    return function() {
      var string0, value1 = value(this), string1;
      if (value1 == null) return void this.removeAttribute(name);
      string0 = this.getAttribute(name);
      string1 = value1 + "";
      return string0 === string1 ? null
          : string0 === string00 && string1 === string10 ? interpolate0
          : (string10 = string1, interpolate0 = interpolate(string00 = string0, value1));
    };
  }

  function attrFunctionNS$1(fullname, interpolate, value) {
    var string00,
        string10,
        interpolate0;
    return function() {
      var string0, value1 = value(this), string1;
      if (value1 == null) return void this.removeAttributeNS(fullname.space, fullname.local);
      string0 = this.getAttributeNS(fullname.space, fullname.local);
      string1 = value1 + "";
      return string0 === string1 ? null
          : string0 === string00 && string1 === string10 ? interpolate0
          : (string10 = string1, interpolate0 = interpolate(string00 = string0, value1));
    };
  }

  function transition_attr(name, value) {
    var fullname = namespace(name), i = fullname === "transform" ? interpolateTransformSvg : interpolate$1;
    return this.attrTween(name, typeof value === "function"
        ? (fullname.local ? attrFunctionNS$1 : attrFunction$1)(fullname, i, tweenValue(this, "attr." + name, value))
        : value == null ? (fullname.local ? attrRemoveNS$1 : attrRemove$1)(fullname)
        : (fullname.local ? attrConstantNS$1 : attrConstant$1)(fullname, i, value));
  }

  function attrInterpolate(name, i) {
    return function(t) {
      this.setAttribute(name, i.call(this, t));
    };
  }

  function attrInterpolateNS(fullname, i) {
    return function(t) {
      this.setAttributeNS(fullname.space, fullname.local, i.call(this, t));
    };
  }

  function attrTweenNS(fullname, value) {
    var t0, i0;
    function tween() {
      var i = value.apply(this, arguments);
      if (i !== i0) t0 = (i0 = i) && attrInterpolateNS(fullname, i);
      return t0;
    }
    tween._value = value;
    return tween;
  }

  function attrTween(name, value) {
    var t0, i0;
    function tween() {
      var i = value.apply(this, arguments);
      if (i !== i0) t0 = (i0 = i) && attrInterpolate(name, i);
      return t0;
    }
    tween._value = value;
    return tween;
  }

  function transition_attrTween(name, value) {
    var key = "attr." + name;
    if (arguments.length < 2) return (key = this.tween(key)) && key._value;
    if (value == null) return this.tween(key, null);
    if (typeof value !== "function") throw new Error;
    var fullname = namespace(name);
    return this.tween(key, (fullname.local ? attrTweenNS : attrTween)(fullname, value));
  }

  function delayFunction(id, value) {
    return function() {
      init(this, id).delay = +value.apply(this, arguments);
    };
  }

  function delayConstant(id, value) {
    return value = +value, function() {
      init(this, id).delay = value;
    };
  }

  function transition_delay(value) {
    var id = this._id;

    return arguments.length
        ? this.each((typeof value === "function"
            ? delayFunction
            : delayConstant)(id, value))
        : get$1(this.node(), id).delay;
  }

  function durationFunction(id, value) {
    return function() {
      set$1(this, id).duration = +value.apply(this, arguments);
    };
  }

  function durationConstant(id, value) {
    return value = +value, function() {
      set$1(this, id).duration = value;
    };
  }

  function transition_duration(value) {
    var id = this._id;

    return arguments.length
        ? this.each((typeof value === "function"
            ? durationFunction
            : durationConstant)(id, value))
        : get$1(this.node(), id).duration;
  }

  function easeConstant(id, value) {
    if (typeof value !== "function") throw new Error;
    return function() {
      set$1(this, id).ease = value;
    };
  }

  function transition_ease(value) {
    var id = this._id;

    return arguments.length
        ? this.each(easeConstant(id, value))
        : get$1(this.node(), id).ease;
  }

  function easeVarying(id, value) {
    return function() {
      var v = value.apply(this, arguments);
      if (typeof v !== "function") throw new Error;
      set$1(this, id).ease = v;
    };
  }

  function transition_easeVarying(value) {
    if (typeof value !== "function") throw new Error;
    return this.each(easeVarying(this._id, value));
  }

  function transition_filter(match) {
    if (typeof match !== "function") match = matcher(match);

    for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, subgroup = subgroups[j] = [], node, i = 0; i < n; ++i) {
        if ((node = group[i]) && match.call(node, node.__data__, i, group)) {
          subgroup.push(node);
        }
      }
    }

    return new Transition(subgroups, this._parents, this._name, this._id);
  }

  function transition_merge(transition) {
    if (transition._id !== this._id) throw new Error;

    for (var groups0 = this._groups, groups1 = transition._groups, m0 = groups0.length, m1 = groups1.length, m = Math.min(m0, m1), merges = new Array(m0), j = 0; j < m; ++j) {
      for (var group0 = groups0[j], group1 = groups1[j], n = group0.length, merge = merges[j] = new Array(n), node, i = 0; i < n; ++i) {
        if (node = group0[i] || group1[i]) {
          merge[i] = node;
        }
      }
    }

    for (; j < m0; ++j) {
      merges[j] = groups0[j];
    }

    return new Transition(merges, this._parents, this._name, this._id);
  }

  function start(name) {
    return (name + "").trim().split(/^|\s+/).every(function(t) {
      var i = t.indexOf(".");
      if (i >= 0) t = t.slice(0, i);
      return !t || t === "start";
    });
  }

  function onFunction(id, name, listener) {
    var on0, on1, sit = start(name) ? init : set$1;
    return function() {
      var schedule = sit(this, id),
          on = schedule.on;

      // If this node shared a dispatch with the previous node,
      // just assign the updated shared dispatch and we’re done!
      // Otherwise, copy-on-write.
      if (on !== on0) (on1 = (on0 = on).copy()).on(name, listener);

      schedule.on = on1;
    };
  }

  function transition_on(name, listener) {
    var id = this._id;

    return arguments.length < 2
        ? get$1(this.node(), id).on.on(name)
        : this.each(onFunction(id, name, listener));
  }

  function removeFunction(id) {
    return function() {
      var parent = this.parentNode;
      for (var i in this.__transition) if (+i !== id) return;
      if (parent) parent.removeChild(this);
    };
  }

  function transition_remove() {
    return this.on("end.remove", removeFunction(this._id));
  }

  function transition_select(select) {
    var name = this._name,
        id = this._id;

    if (typeof select !== "function") select = selector(select);

    for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, subgroup = subgroups[j] = new Array(n), node, subnode, i = 0; i < n; ++i) {
        if ((node = group[i]) && (subnode = select.call(node, node.__data__, i, group))) {
          if ("__data__" in node) subnode.__data__ = node.__data__;
          subgroup[i] = subnode;
          schedule(subgroup[i], name, id, i, subgroup, get$1(node, id));
        }
      }
    }

    return new Transition(subgroups, this._parents, name, id);
  }

  function transition_selectAll(select) {
    var name = this._name,
        id = this._id;

    if (typeof select !== "function") select = selectorAll(select);

    for (var groups = this._groups, m = groups.length, subgroups = [], parents = [], j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
        if (node = group[i]) {
          for (var children = select.call(node, node.__data__, i, group), child, inherit = get$1(node, id), k = 0, l = children.length; k < l; ++k) {
            if (child = children[k]) {
              schedule(child, name, id, k, children, inherit);
            }
          }
          subgroups.push(children);
          parents.push(node);
        }
      }
    }

    return new Transition(subgroups, parents, name, id);
  }

  var Selection$1 = selection.prototype.constructor;

  function transition_selection() {
    return new Selection$1(this._groups, this._parents);
  }

  function styleNull(name, interpolate) {
    var string00,
        string10,
        interpolate0;
    return function() {
      var string0 = styleValue(this, name),
          string1 = (this.style.removeProperty(name), styleValue(this, name));
      return string0 === string1 ? null
          : string0 === string00 && string1 === string10 ? interpolate0
          : interpolate0 = interpolate(string00 = string0, string10 = string1);
    };
  }

  function styleRemove$1(name) {
    return function() {
      this.style.removeProperty(name);
    };
  }

  function styleConstant$1(name, interpolate, value1) {
    var string00,
        string1 = value1 + "",
        interpolate0;
    return function() {
      var string0 = styleValue(this, name);
      return string0 === string1 ? null
          : string0 === string00 ? interpolate0
          : interpolate0 = interpolate(string00 = string0, value1);
    };
  }

  function styleFunction$1(name, interpolate, value) {
    var string00,
        string10,
        interpolate0;
    return function() {
      var string0 = styleValue(this, name),
          value1 = value(this),
          string1 = value1 + "";
      if (value1 == null) string1 = value1 = (this.style.removeProperty(name), styleValue(this, name));
      return string0 === string1 ? null
          : string0 === string00 && string1 === string10 ? interpolate0
          : (string10 = string1, interpolate0 = interpolate(string00 = string0, value1));
    };
  }

  function styleMaybeRemove(id, name) {
    var on0, on1, listener0, key = "style." + name, event = "end." + key, remove;
    return function() {
      var schedule = set$1(this, id),
          on = schedule.on,
          listener = schedule.value[key] == null ? remove || (remove = styleRemove$1(name)) : undefined;

      // If this node shared a dispatch with the previous node,
      // just assign the updated shared dispatch and we’re done!
      // Otherwise, copy-on-write.
      if (on !== on0 || listener0 !== listener) (on1 = (on0 = on).copy()).on(event, listener0 = listener);

      schedule.on = on1;
    };
  }

  function transition_style(name, value, priority) {
    var i = (name += "") === "transform" ? interpolateTransformCss : interpolate$1;
    return value == null ? this
        .styleTween(name, styleNull(name, i))
        .on("end.style." + name, styleRemove$1(name))
      : typeof value === "function" ? this
        .styleTween(name, styleFunction$1(name, i, tweenValue(this, "style." + name, value)))
        .each(styleMaybeRemove(this._id, name))
      : this
        .styleTween(name, styleConstant$1(name, i, value), priority)
        .on("end.style." + name, null);
  }

  function styleInterpolate(name, i, priority) {
    return function(t) {
      this.style.setProperty(name, i.call(this, t), priority);
    };
  }

  function styleTween(name, value, priority) {
    var t, i0;
    function tween() {
      var i = value.apply(this, arguments);
      if (i !== i0) t = (i0 = i) && styleInterpolate(name, i, priority);
      return t;
    }
    tween._value = value;
    return tween;
  }

  function transition_styleTween(name, value, priority) {
    var key = "style." + (name += "");
    if (arguments.length < 2) return (key = this.tween(key)) && key._value;
    if (value == null) return this.tween(key, null);
    if (typeof value !== "function") throw new Error;
    return this.tween(key, styleTween(name, value, priority == null ? "" : priority));
  }

  function textConstant$1(value) {
    return function() {
      this.textContent = value;
    };
  }

  function textFunction$1(value) {
    return function() {
      var value1 = value(this);
      this.textContent = value1 == null ? "" : value1;
    };
  }

  function transition_text(value) {
    return this.tween("text", typeof value === "function"
        ? textFunction$1(tweenValue(this, "text", value))
        : textConstant$1(value == null ? "" : value + ""));
  }

  function textInterpolate(i) {
    return function(t) {
      this.textContent = i.call(this, t);
    };
  }

  function textTween(value) {
    var t0, i0;
    function tween() {
      var i = value.apply(this, arguments);
      if (i !== i0) t0 = (i0 = i) && textInterpolate(i);
      return t0;
    }
    tween._value = value;
    return tween;
  }

  function transition_textTween(value) {
    var key = "text";
    if (arguments.length < 1) return (key = this.tween(key)) && key._value;
    if (value == null) return this.tween(key, null);
    if (typeof value !== "function") throw new Error;
    return this.tween(key, textTween(value));
  }

  function transition_transition() {
    var name = this._name,
        id0 = this._id,
        id1 = newId();

    for (var groups = this._groups, m = groups.length, j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
        if (node = group[i]) {
          var inherit = get$1(node, id0);
          schedule(node, name, id1, i, group, {
            time: inherit.time + inherit.delay + inherit.duration,
            delay: 0,
            duration: inherit.duration,
            ease: inherit.ease
          });
        }
      }
    }

    return new Transition(groups, this._parents, name, id1);
  }

  function transition_end() {
    var on0, on1, that = this, id = that._id, size = that.size();
    return new Promise(function(resolve, reject) {
      var cancel = {value: reject},
          end = {value: function() { if (--size === 0) resolve(); }};

      that.each(function() {
        var schedule = set$1(this, id),
            on = schedule.on;

        // If this node shared a dispatch with the previous node,
        // just assign the updated shared dispatch and we’re done!
        // Otherwise, copy-on-write.
        if (on !== on0) {
          on1 = (on0 = on).copy();
          on1._.cancel.push(cancel);
          on1._.interrupt.push(cancel);
          on1._.end.push(end);
        }

        schedule.on = on1;
      });

      // The selection was empty, resolve end immediately
      if (size === 0) resolve();
    });
  }

  var id = 0;

  function Transition(groups, parents, name, id) {
    this._groups = groups;
    this._parents = parents;
    this._name = name;
    this._id = id;
  }

  function transition(name) {
    return selection().transition(name);
  }

  function newId() {
    return ++id;
  }

  var selection_prototype = selection.prototype;

  Transition.prototype = transition.prototype = {
    constructor: Transition,
    select: transition_select,
    selectAll: transition_selectAll,
    selectChild: selection_prototype.selectChild,
    selectChildren: selection_prototype.selectChildren,
    filter: transition_filter,
    merge: transition_merge,
    selection: transition_selection,
    transition: transition_transition,
    call: selection_prototype.call,
    nodes: selection_prototype.nodes,
    node: selection_prototype.node,
    size: selection_prototype.size,
    empty: selection_prototype.empty,
    each: selection_prototype.each,
    on: transition_on,
    attr: transition_attr,
    attrTween: transition_attrTween,
    style: transition_style,
    styleTween: transition_styleTween,
    text: transition_text,
    textTween: transition_textTween,
    remove: transition_remove,
    tween: transition_tween,
    delay: transition_delay,
    duration: transition_duration,
    ease: transition_ease,
    easeVarying: transition_easeVarying,
    end: transition_end,
    [Symbol.iterator]: selection_prototype[Symbol.iterator]
  };

  function cubicInOut(t) {
    return ((t *= 2) <= 1 ? t * t * t : (t -= 2) * t * t + 2) / 2;
  }

  var defaultTiming = {
    time: null, // Set on use.
    delay: 0,
    duration: 250,
    ease: cubicInOut
  };

  function inherit(node, id) {
    var timing;
    while (!(timing = node.__transition) || !(timing = timing[id])) {
      if (!(node = node.parentNode)) {
        throw new Error(`transition ${id} not found`);
      }
    }
    return timing;
  }

  function selection_transition(name) {
    var id,
        timing;

    if (name instanceof Transition) {
      id = name._id, name = name._name;
    } else {
      id = newId(), (timing = defaultTiming).time = now(), name = name == null ? null : name + "";
    }

    for (var groups = this._groups, m = groups.length, j = 0; j < m; ++j) {
      for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
        if (node = group[i]) {
          schedule(node, name, id, i, group, timing || inherit(node, id));
        }
      }
    }

    return new Transition(groups, this._parents, name, id);
  }

  selection.prototype.interrupt = selection_interrupt;
  selection.prototype.transition = selection_transition;

  var constant$3 = x => () => x;

  function ZoomEvent(type, {
    sourceEvent,
    target,
    transform,
    dispatch
  }) {
    Object.defineProperties(this, {
      type: {value: type, enumerable: true, configurable: true},
      sourceEvent: {value: sourceEvent, enumerable: true, configurable: true},
      target: {value: target, enumerable: true, configurable: true},
      transform: {value: transform, enumerable: true, configurable: true},
      _: {value: dispatch}
    });
  }

  function Transform(k, x, y) {
    this.k = k;
    this.x = x;
    this.y = y;
  }

  Transform.prototype = {
    constructor: Transform,
    scale: function(k) {
      return k === 1 ? this : new Transform(this.k * k, this.x, this.y);
    },
    translate: function(x, y) {
      return x === 0 & y === 0 ? this : new Transform(this.k, this.x + this.k * x, this.y + this.k * y);
    },
    apply: function(point) {
      return [point[0] * this.k + this.x, point[1] * this.k + this.y];
    },
    applyX: function(x) {
      return x * this.k + this.x;
    },
    applyY: function(y) {
      return y * this.k + this.y;
    },
    invert: function(location) {
      return [(location[0] - this.x) / this.k, (location[1] - this.y) / this.k];
    },
    invertX: function(x) {
      return (x - this.x) / this.k;
    },
    invertY: function(y) {
      return (y - this.y) / this.k;
    },
    rescaleX: function(x) {
      return x.copy().domain(x.range().map(this.invertX, this).map(x.invert, x));
    },
    rescaleY: function(y) {
      return y.copy().domain(y.range().map(this.invertY, this).map(y.invert, y));
    },
    toString: function() {
      return "translate(" + this.x + "," + this.y + ") scale(" + this.k + ")";
    }
  };

  var identity$4 = new Transform(1, 0, 0);

  function nopropagation(event) {
    event.stopImmediatePropagation();
  }

  function noevent$1(event) {
    event.preventDefault();
    event.stopImmediatePropagation();
  }

  // Ignore right-click, since that should open the context menu.
  // except for pinch-to-zoom, which is sent as a wheel+ctrlKey event
  function defaultFilter(event) {
    return (!event.ctrlKey || event.type === 'wheel') && !event.button;
  }

  function defaultExtent() {
    var e = this;
    if (e instanceof SVGElement) {
      e = e.ownerSVGElement || e;
      if (e.hasAttribute("viewBox")) {
        e = e.viewBox.baseVal;
        return [[e.x, e.y], [e.x + e.width, e.y + e.height]];
      }
      return [[0, 0], [e.width.baseVal.value, e.height.baseVal.value]];
    }
    return [[0, 0], [e.clientWidth, e.clientHeight]];
  }

  function defaultTransform() {
    return this.__zoom || identity$4;
  }

  function defaultWheelDelta(event) {
    return -event.deltaY * (event.deltaMode === 1 ? 0.05 : event.deltaMode ? 1 : 0.002) * (event.ctrlKey ? 10 : 1);
  }

  function defaultTouchable() {
    return navigator.maxTouchPoints || ("ontouchstart" in this);
  }

  function defaultConstrain(transform, extent, translateExtent) {
    var dx0 = transform.invertX(extent[0][0]) - translateExtent[0][0],
        dx1 = transform.invertX(extent[1][0]) - translateExtent[1][0],
        dy0 = transform.invertY(extent[0][1]) - translateExtent[0][1],
        dy1 = transform.invertY(extent[1][1]) - translateExtent[1][1];
    return transform.translate(
      dx1 > dx0 ? (dx0 + dx1) / 2 : Math.min(0, dx0) || Math.max(0, dx1),
      dy1 > dy0 ? (dy0 + dy1) / 2 : Math.min(0, dy0) || Math.max(0, dy1)
    );
  }

  function zoom() {
    var filter = defaultFilter,
        extent = defaultExtent,
        constrain = defaultConstrain,
        wheelDelta = defaultWheelDelta,
        touchable = defaultTouchable,
        scaleExtent = [0, Infinity],
        translateExtent = [[-Infinity, -Infinity], [Infinity, Infinity]],
        duration = 250,
        interpolate = interpolateZoom,
        listeners = dispatch("start", "zoom", "end"),
        touchstarting,
        touchfirst,
        touchending,
        touchDelay = 500,
        wheelDelay = 150,
        clickDistance2 = 0,
        tapDistance = 10;

    function zoom(selection) {
      selection
          .property("__zoom", defaultTransform)
          .on("wheel.zoom", wheeled, {passive: false})
          .on("mousedown.zoom", mousedowned)
          .on("dblclick.zoom", dblclicked)
        .filter(touchable)
          .on("touchstart.zoom", touchstarted)
          .on("touchmove.zoom", touchmoved)
          .on("touchend.zoom touchcancel.zoom", touchended)
          .style("-webkit-tap-highlight-color", "rgba(0,0,0,0)");
    }

    zoom.transform = function(collection, transform, point, event) {
      var selection = collection.selection ? collection.selection() : collection;
      selection.property("__zoom", defaultTransform);
      if (collection !== selection) {
        schedule(collection, transform, point, event);
      } else {
        selection.interrupt().each(function() {
          gesture(this, arguments)
            .event(event)
            .start()
            .zoom(null, typeof transform === "function" ? transform.apply(this, arguments) : transform)
            .end();
        });
      }
    };

    zoom.scaleBy = function(selection, k, p, event) {
      zoom.scaleTo(selection, function() {
        var k0 = this.__zoom.k,
            k1 = typeof k === "function" ? k.apply(this, arguments) : k;
        return k0 * k1;
      }, p, event);
    };

    zoom.scaleTo = function(selection, k, p, event) {
      zoom.transform(selection, function() {
        var e = extent.apply(this, arguments),
            t0 = this.__zoom,
            p0 = p == null ? centroid(e) : typeof p === "function" ? p.apply(this, arguments) : p,
            p1 = t0.invert(p0),
            k1 = typeof k === "function" ? k.apply(this, arguments) : k;
        return constrain(translate(scale(t0, k1), p0, p1), e, translateExtent);
      }, p, event);
    };

    zoom.translateBy = function(selection, x, y, event) {
      zoom.transform(selection, function() {
        return constrain(this.__zoom.translate(
          typeof x === "function" ? x.apply(this, arguments) : x,
          typeof y === "function" ? y.apply(this, arguments) : y
        ), extent.apply(this, arguments), translateExtent);
      }, null, event);
    };

    zoom.translateTo = function(selection, x, y, p, event) {
      zoom.transform(selection, function() {
        var e = extent.apply(this, arguments),
            t = this.__zoom,
            p0 = p == null ? centroid(e) : typeof p === "function" ? p.apply(this, arguments) : p;
        return constrain(identity$4.translate(p0[0], p0[1]).scale(t.k).translate(
          typeof x === "function" ? -x.apply(this, arguments) : -x,
          typeof y === "function" ? -y.apply(this, arguments) : -y
        ), e, translateExtent);
      }, p, event);
    };

    function scale(transform, k) {
      k = Math.max(scaleExtent[0], Math.min(scaleExtent[1], k));
      return k === transform.k ? transform : new Transform(k, transform.x, transform.y);
    }

    function translate(transform, p0, p1) {
      var x = p0[0] - p1[0] * transform.k, y = p0[1] - p1[1] * transform.k;
      return x === transform.x && y === transform.y ? transform : new Transform(transform.k, x, y);
    }

    function centroid(extent) {
      return [(+extent[0][0] + +extent[1][0]) / 2, (+extent[0][1] + +extent[1][1]) / 2];
    }

    function schedule(transition, transform, point, event) {
      transition
          .on("start.zoom", function() { gesture(this, arguments).event(event).start(); })
          .on("interrupt.zoom end.zoom", function() { gesture(this, arguments).event(event).end(); })
          .tween("zoom", function() {
            var that = this,
                args = arguments,
                g = gesture(that, args).event(event),
                e = extent.apply(that, args),
                p = point == null ? centroid(e) : typeof point === "function" ? point.apply(that, args) : point,
                w = Math.max(e[1][0] - e[0][0], e[1][1] - e[0][1]),
                a = that.__zoom,
                b = typeof transform === "function" ? transform.apply(that, args) : transform,
                i = interpolate(a.invert(p).concat(w / a.k), b.invert(p).concat(w / b.k));
            return function(t) {
              if (t === 1) t = b; // Avoid rounding error on end.
              else { var l = i(t), k = w / l[2]; t = new Transform(k, p[0] - l[0] * k, p[1] - l[1] * k); }
              g.zoom(null, t);
            };
          });
    }

    function gesture(that, args, clean) {
      return (!clean && that.__zooming) || new Gesture(that, args);
    }

    function Gesture(that, args) {
      this.that = that;
      this.args = args;
      this.active = 0;
      this.sourceEvent = null;
      this.extent = extent.apply(that, args);
      this.taps = 0;
    }

    Gesture.prototype = {
      event: function(event) {
        if (event) this.sourceEvent = event;
        return this;
      },
      start: function() {
        if (++this.active === 1) {
          this.that.__zooming = this;
          this.emit("start");
        }
        return this;
      },
      zoom: function(key, transform) {
        if (this.mouse && key !== "mouse") this.mouse[1] = transform.invert(this.mouse[0]);
        if (this.touch0 && key !== "touch") this.touch0[1] = transform.invert(this.touch0[0]);
        if (this.touch1 && key !== "touch") this.touch1[1] = transform.invert(this.touch1[0]);
        this.that.__zoom = transform;
        this.emit("zoom");
        return this;
      },
      end: function() {
        if (--this.active === 0) {
          delete this.that.__zooming;
          this.emit("end");
        }
        return this;
      },
      emit: function(type) {
        var d = select(this.that).datum();
        listeners.call(
          type,
          this.that,
          new ZoomEvent(type, {
            sourceEvent: this.sourceEvent,
            target: zoom,
            type,
            transform: this.that.__zoom,
            dispatch: listeners
          }),
          d
        );
      }
    };

    function wheeled(event, ...args) {
      if (!filter.apply(this, arguments)) return;
      var g = gesture(this, args).event(event),
          t = this.__zoom,
          k = Math.max(scaleExtent[0], Math.min(scaleExtent[1], t.k * Math.pow(2, wheelDelta.apply(this, arguments)))),
          p = pointer(event);

      // If the mouse is in the same location as before, reuse it.
      // If there were recent wheel events, reset the wheel idle timeout.
      if (g.wheel) {
        if (g.mouse[0][0] !== p[0] || g.mouse[0][1] !== p[1]) {
          g.mouse[1] = t.invert(g.mouse[0] = p);
        }
        clearTimeout(g.wheel);
      }

      // If this wheel event won’t trigger a transform change, ignore it.
      else if (t.k === k) return;

      // Otherwise, capture the mouse point and location at the start.
      else {
        g.mouse = [p, t.invert(p)];
        interrupt(this);
        g.start();
      }

      noevent$1(event);
      g.wheel = setTimeout(wheelidled, wheelDelay);
      g.zoom("mouse", constrain(translate(scale(t, k), g.mouse[0], g.mouse[1]), g.extent, translateExtent));

      function wheelidled() {
        g.wheel = null;
        g.end();
      }
    }

    function mousedowned(event, ...args) {
      if (touchending || !filter.apply(this, arguments)) return;
      var currentTarget = event.currentTarget,
          g = gesture(this, args, true).event(event),
          v = select(event.view).on("mousemove.zoom", mousemoved, true).on("mouseup.zoom", mouseupped, true),
          p = pointer(event, currentTarget),
          x0 = event.clientX,
          y0 = event.clientY;

      dragDisable(event.view);
      nopropagation(event);
      g.mouse = [p, this.__zoom.invert(p)];
      interrupt(this);
      g.start();

      function mousemoved(event) {
        noevent$1(event);
        if (!g.moved) {
          var dx = event.clientX - x0, dy = event.clientY - y0;
          g.moved = dx * dx + dy * dy > clickDistance2;
        }
        g.event(event)
         .zoom("mouse", constrain(translate(g.that.__zoom, g.mouse[0] = pointer(event, currentTarget), g.mouse[1]), g.extent, translateExtent));
      }

      function mouseupped(event) {
        v.on("mousemove.zoom mouseup.zoom", null);
        yesdrag(event.view, g.moved);
        noevent$1(event);
        g.event(event).end();
      }
    }

    function dblclicked(event, ...args) {
      if (!filter.apply(this, arguments)) return;
      var t0 = this.__zoom,
          p0 = pointer(event.changedTouches ? event.changedTouches[0] : event, this),
          p1 = t0.invert(p0),
          k1 = t0.k * (event.shiftKey ? 0.5 : 2),
          t1 = constrain(translate(scale(t0, k1), p0, p1), extent.apply(this, args), translateExtent);

      noevent$1(event);
      if (duration > 0) select(this).transition().duration(duration).call(schedule, t1, p0, event);
      else select(this).call(zoom.transform, t1, p0, event);
    }

    function touchstarted(event, ...args) {
      if (!filter.apply(this, arguments)) return;
      var touches = event.touches,
          n = touches.length,
          g = gesture(this, args, event.changedTouches.length === n).event(event),
          started, i, t, p;

      nopropagation(event);
      for (i = 0; i < n; ++i) {
        t = touches[i], p = pointer(t, this);
        p = [p, this.__zoom.invert(p), t.identifier];
        if (!g.touch0) g.touch0 = p, started = true, g.taps = 1 + !!touchstarting;
        else if (!g.touch1 && g.touch0[2] !== p[2]) g.touch1 = p, g.taps = 0;
      }

      if (touchstarting) touchstarting = clearTimeout(touchstarting);

      if (started) {
        if (g.taps < 2) touchfirst = p[0], touchstarting = setTimeout(function() { touchstarting = null; }, touchDelay);
        interrupt(this);
        g.start();
      }
    }

    function touchmoved(event, ...args) {
      if (!this.__zooming) return;
      var g = gesture(this, args).event(event),
          touches = event.changedTouches,
          n = touches.length, i, t, p, l;

      noevent$1(event);
      for (i = 0; i < n; ++i) {
        t = touches[i], p = pointer(t, this);
        if (g.touch0 && g.touch0[2] === t.identifier) g.touch0[0] = p;
        else if (g.touch1 && g.touch1[2] === t.identifier) g.touch1[0] = p;
      }
      t = g.that.__zoom;
      if (g.touch1) {
        var p0 = g.touch0[0], l0 = g.touch0[1],
            p1 = g.touch1[0], l1 = g.touch1[1],
            dp = (dp = p1[0] - p0[0]) * dp + (dp = p1[1] - p0[1]) * dp,
            dl = (dl = l1[0] - l0[0]) * dl + (dl = l1[1] - l0[1]) * dl;
        t = scale(t, Math.sqrt(dp / dl));
        p = [(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2];
        l = [(l0[0] + l1[0]) / 2, (l0[1] + l1[1]) / 2];
      }
      else if (g.touch0) p = g.touch0[0], l = g.touch0[1];
      else return;

      g.zoom("touch", constrain(translate(t, p, l), g.extent, translateExtent));
    }

    function touchended(event, ...args) {
      if (!this.__zooming) return;
      var g = gesture(this, args).event(event),
          touches = event.changedTouches,
          n = touches.length, i, t;

      nopropagation(event);
      if (touchending) clearTimeout(touchending);
      touchending = setTimeout(function() { touchending = null; }, touchDelay);
      for (i = 0; i < n; ++i) {
        t = touches[i];
        if (g.touch0 && g.touch0[2] === t.identifier) delete g.touch0;
        else if (g.touch1 && g.touch1[2] === t.identifier) delete g.touch1;
      }
      if (g.touch1 && !g.touch0) g.touch0 = g.touch1, delete g.touch1;
      if (g.touch0) g.touch0[1] = this.__zoom.invert(g.touch0[0]);
      else {
        g.end();
        // If this was a dbltap, reroute to the (optional) dblclick.zoom handler.
        if (g.taps === 2) {
          t = pointer(t, this);
          if (Math.hypot(touchfirst[0] - t[0], touchfirst[1] - t[1]) < tapDistance) {
            var p = select(this).on("dblclick.zoom");
            if (p) p.apply(this, arguments);
          }
        }
      }
    }

    zoom.wheelDelta = function(_) {
      return arguments.length ? (wheelDelta = typeof _ === "function" ? _ : constant$3(+_), zoom) : wheelDelta;
    };

    zoom.filter = function(_) {
      return arguments.length ? (filter = typeof _ === "function" ? _ : constant$3(!!_), zoom) : filter;
    };

    zoom.touchable = function(_) {
      return arguments.length ? (touchable = typeof _ === "function" ? _ : constant$3(!!_), zoom) : touchable;
    };

    zoom.extent = function(_) {
      return arguments.length ? (extent = typeof _ === "function" ? _ : constant$3([[+_[0][0], +_[0][1]], [+_[1][0], +_[1][1]]]), zoom) : extent;
    };

    zoom.scaleExtent = function(_) {
      return arguments.length ? (scaleExtent[0] = +_[0], scaleExtent[1] = +_[1], zoom) : [scaleExtent[0], scaleExtent[1]];
    };

    zoom.translateExtent = function(_) {
      return arguments.length ? (translateExtent[0][0] = +_[0][0], translateExtent[1][0] = +_[1][0], translateExtent[0][1] = +_[0][1], translateExtent[1][1] = +_[1][1], zoom) : [[translateExtent[0][0], translateExtent[0][1]], [translateExtent[1][0], translateExtent[1][1]]];
    };

    zoom.constrain = function(_) {
      return arguments.length ? (constrain = _, zoom) : constrain;
    };

    zoom.duration = function(_) {
      return arguments.length ? (duration = +_, zoom) : duration;
    };

    zoom.interpolate = function(_) {
      return arguments.length ? (interpolate = _, zoom) : interpolate;
    };

    zoom.on = function() {
      var value = listeners.on.apply(listeners, arguments);
      return value === listeners ? zoom : value;
    };

    zoom.clickDistance = function(_) {
      return arguments.length ? (clickDistance2 = (_ = +_) * _, zoom) : Math.sqrt(clickDistance2);
    };

    zoom.tapDistance = function(_) {
      return arguments.length ? (tapDistance = +_, zoom) : tapDistance;
    };

    return zoom;
  }

  function defaultSeparation(a, b) {
    return a.parent === b.parent ? 1 : 2;
  }

  function meanX(children) {
    return children.reduce(meanXReduce, 0) / children.length;
  }

  function meanXReduce(x, c) {
    return x + c.x;
  }

  function maxY(children) {
    return 1 + children.reduce(maxYReduce, 0);
  }

  function maxYReduce(y, c) {
    return Math.max(y, c.y);
  }

  function leafLeft(node) {
    var children;
    while (children = node.children) node = children[0];
    return node;
  }

  function leafRight(node) {
    var children;
    while (children = node.children) node = children[children.length - 1];
    return node;
  }

  function cluster() {
    var separation = defaultSeparation,
        dx = 1,
        dy = 1,
        nodeSize = false;

    function cluster(root) {
      var previousNode,
          x = 0;

      // First walk, computing the initial x & y values.
      root.eachAfter(function(node) {
        var children = node.children;
        if (children) {
          node.x = meanX(children);
          node.y = maxY(children);
        } else {
          node.x = previousNode ? x += separation(node, previousNode) : 0;
          node.y = 0;
          previousNode = node;
        }
      });

      var left = leafLeft(root),
          right = leafRight(root),
          x0 = left.x - separation(left, right) / 2,
          x1 = right.x + separation(right, left) / 2;

      // Second walk, normalizing x & y to the desired size.
      return root.eachAfter(nodeSize ? function(node) {
        node.x = (node.x - root.x) * dx;
        node.y = (root.y - node.y) * dy;
      } : function(node) {
        node.x = (node.x - x0) / (x1 - x0) * dx;
        node.y = (1 - (root.y ? node.y / root.y : 1)) * dy;
      });
    }

    cluster.separation = function(x) {
      return arguments.length ? (separation = x, cluster) : separation;
    };

    cluster.size = function(x) {
      return arguments.length ? (nodeSize = false, dx = +x[0], dy = +x[1], cluster) : (nodeSize ? null : [dx, dy]);
    };

    cluster.nodeSize = function(x) {
      return arguments.length ? (nodeSize = true, dx = +x[0], dy = +x[1], cluster) : (nodeSize ? [dx, dy] : null);
    };

    return cluster;
  }

  function count(node) {
    var sum = 0,
        children = node.children,
        i = children && children.length;
    if (!i) sum = 1;
    else while (--i >= 0) sum += children[i].value;
    node.value = sum;
  }

  function node_count() {
    return this.eachAfter(count);
  }

  function node_each(callback, that) {
    let index = -1;
    for (const node of this) {
      callback.call(that, node, ++index, this);
    }
    return this;
  }

  function node_eachBefore(callback, that) {
    var node = this, nodes = [node], children, i, index = -1;
    while (node = nodes.pop()) {
      callback.call(that, node, ++index, this);
      if (children = node.children) {
        for (i = children.length - 1; i >= 0; --i) {
          nodes.push(children[i]);
        }
      }
    }
    return this;
  }

  function node_eachAfter(callback, that) {
    var node = this, nodes = [node], next = [], children, i, n, index = -1;
    while (node = nodes.pop()) {
      next.push(node);
      if (children = node.children) {
        for (i = 0, n = children.length; i < n; ++i) {
          nodes.push(children[i]);
        }
      }
    }
    while (node = next.pop()) {
      callback.call(that, node, ++index, this);
    }
    return this;
  }

  function node_find(callback, that) {
    let index = -1;
    for (const node of this) {
      if (callback.call(that, node, ++index, this)) {
        return node;
      }
    }
  }

  function node_sum(value) {
    return this.eachAfter(function(node) {
      var sum = +value(node.data) || 0,
          children = node.children,
          i = children && children.length;
      while (--i >= 0) sum += children[i].value;
      node.value = sum;
    });
  }

  function node_sort(compare) {
    return this.eachBefore(function(node) {
      if (node.children) {
        node.children.sort(compare);
      }
    });
  }

  function node_path(end) {
    var start = this,
        ancestor = leastCommonAncestor(start, end),
        nodes = [start];
    while (start !== ancestor) {
      start = start.parent;
      nodes.push(start);
    }
    var k = nodes.length;
    while (end !== ancestor) {
      nodes.splice(k, 0, end);
      end = end.parent;
    }
    return nodes;
  }

  function leastCommonAncestor(a, b) {
    if (a === b) return a;
    var aNodes = a.ancestors(),
        bNodes = b.ancestors(),
        c = null;
    a = aNodes.pop();
    b = bNodes.pop();
    while (a === b) {
      c = a;
      a = aNodes.pop();
      b = bNodes.pop();
    }
    return c;
  }

  function node_ancestors() {
    var node = this, nodes = [node];
    while (node = node.parent) {
      nodes.push(node);
    }
    return nodes;
  }

  function node_descendants() {
    return Array.from(this);
  }

  function node_leaves() {
    var leaves = [];
    this.eachBefore(function(node) {
      if (!node.children) {
        leaves.push(node);
      }
    });
    return leaves;
  }

  function node_links() {
    var root = this, links = [];
    root.each(function(node) {
      if (node !== root) { // Don’t include the root’s parent, if any.
        links.push({source: node.parent, target: node});
      }
    });
    return links;
  }

  function* node_iterator() {
    var node = this, current, next = [node], children, i, n;
    do {
      current = next.reverse(), next = [];
      while (node = current.pop()) {
        yield node;
        if (children = node.children) {
          for (i = 0, n = children.length; i < n; ++i) {
            next.push(children[i]);
          }
        }
      }
    } while (next.length);
  }

  function hierarchy(data, children) {
    if (data instanceof Map) {
      data = [undefined, data];
      if (children === undefined) children = mapChildren;
    } else if (children === undefined) {
      children = objectChildren;
    }

    var root = new Node(data),
        node,
        nodes = [root],
        child,
        childs,
        i,
        n;

    while (node = nodes.pop()) {
      if ((childs = children(node.data)) && (n = (childs = Array.from(childs)).length)) {
        node.children = childs;
        for (i = n - 1; i >= 0; --i) {
          nodes.push(child = childs[i] = new Node(childs[i]));
          child.parent = node;
          child.depth = node.depth + 1;
        }
      }
    }

    return root.eachBefore(computeHeight);
  }

  function node_copy() {
    return hierarchy(this).eachBefore(copyData);
  }

  function objectChildren(d) {
    return d.children;
  }

  function mapChildren(d) {
    return Array.isArray(d) ? d[1] : null;
  }

  function copyData(node) {
    if (node.data.value !== undefined) node.value = node.data.value;
    node.data = node.data.data;
  }

  function computeHeight(node) {
    var height = 0;
    do node.height = height;
    while ((node = node.parent) && (node.height < ++height));
  }

  function Node(data) {
    this.data = data;
    this.depth =
    this.height = 0;
    this.parent = null;
  }

  Node.prototype = hierarchy.prototype = {
    constructor: Node,
    count: node_count,
    each: node_each,
    eachAfter: node_eachAfter,
    eachBefore: node_eachBefore,
    find: node_find,
    sum: node_sum,
    sort: node_sort,
    path: node_path,
    ancestors: node_ancestors,
    descendants: node_descendants,
    leaves: node_leaves,
    links: node_links,
    copy: node_copy,
    [Symbol.iterator]: node_iterator
  };

  function defaultSeparation$1(a, b) {
    return a.parent === b.parent ? 1 : 2;
  }

  // function radialSeparation(a, b) {
  //   return (a.parent === b.parent ? 1 : 2) / a.depth;
  // }

  // This function is used to traverse the left contour of a subtree (or
  // subforest). It returns the successor of v on this contour. This successor is
  // either given by the leftmost child of v or by the thread of v. The function
  // returns null if and only if v is on the highest level of its subtree.
  function nextLeft(v) {
    var children = v.children;
    return children ? children[0] : v.t;
  }

  // This function works analogously to nextLeft.
  function nextRight(v) {
    var children = v.children;
    return children ? children[children.length - 1] : v.t;
  }

  // Shifts the current subtree rooted at w+. This is done by increasing
  // prelim(w+) and mod(w+) by shift.
  function moveSubtree(wm, wp, shift) {
    var change = shift / (wp.i - wm.i);
    wp.c -= change;
    wp.s += shift;
    wm.c += change;
    wp.z += shift;
    wp.m += shift;
  }

  // All other shifts, applied to the smaller subtrees between w- and w+, are
  // performed by this function. To prepare the shifts, we have to adjust
  // change(w+), shift(w+), and change(w-).
  function executeShifts(v) {
    var shift = 0,
        change = 0,
        children = v.children,
        i = children.length,
        w;
    while (--i >= 0) {
      w = children[i];
      w.z += shift;
      w.m += shift;
      shift += w.s + (change += w.c);
    }
  }

  // If vi-’s ancestor is a sibling of v, returns vi-’s ancestor. Otherwise,
  // returns the specified (default) ancestor.
  function nextAncestor(vim, v, ancestor) {
    return vim.a.parent === v.parent ? vim.a : ancestor;
  }

  function TreeNode(node, i) {
    this._ = node;
    this.parent = null;
    this.children = null;
    this.A = null; // default ancestor
    this.a = this; // ancestor
    this.z = 0; // prelim
    this.m = 0; // mod
    this.c = 0; // change
    this.s = 0; // shift
    this.t = null; // thread
    this.i = i; // number
  }

  TreeNode.prototype = Object.create(Node.prototype);

  function treeRoot(root) {
    var tree = new TreeNode(root, 0),
        node,
        nodes = [tree],
        child,
        children,
        i,
        n;

    while (node = nodes.pop()) {
      if (children = node._.children) {
        node.children = new Array(n = children.length);
        for (i = n - 1; i >= 0; --i) {
          nodes.push(child = node.children[i] = new TreeNode(children[i], i));
          child.parent = node;
        }
      }
    }

    (tree.parent = new TreeNode(null, 0)).children = [tree];
    return tree;
  }

  // Node-link tree diagram using the Reingold-Tilford "tidy" algorithm
  function d3Tree() {
    var separation = defaultSeparation$1,
        dx = 1,
        dy = 1,
        nodeSize = null;

    function tree(root) {
      var t = treeRoot(root);

      // Compute the layout using Buchheim et al.’s algorithm.
      t.eachAfter(firstWalk), t.parent.m = -t.z;
      t.eachBefore(secondWalk);

      // If a fixed node size is specified, scale x and y.
      if (nodeSize) root.eachBefore(sizeNode);

      // If a fixed tree size is specified, scale x and y based on the extent.
      // Compute the left-most, right-most, and depth-most nodes for extents.
      else {
        var left = root,
            right = root,
            bottom = root;
        root.eachBefore(function(node) {
          if (node.x < left.x) left = node;
          if (node.x > right.x) right = node;
          if (node.depth > bottom.depth) bottom = node;
        });
        var s = left === right ? 1 : separation(left, right) / 2,
            tx = s - left.x,
            kx = dx / (right.x + s + tx),
            ky = dy / (bottom.depth || 1);
        root.eachBefore(function(node) {
          node.x = (node.x + tx) * kx;
          node.y = node.depth * ky;
        });
      }

      return root;
    }

    // Computes a preliminary x-coordinate for v. Before that, FIRST WALK is
    // applied recursively to the children of v, as well as the function
    // APPORTION. After spacing out the children by calling EXECUTE SHIFTS, the
    // node v is placed to the midpoint of its outermost children.
    function firstWalk(v) {
      var children = v.children,
          siblings = v.parent.children,
          w = v.i ? siblings[v.i - 1] : null;
      if (children) {
        executeShifts(v);
        var midpoint = (children[0].z + children[children.length - 1].z) / 2;
        if (w) {
          v.z = w.z + separation(v._, w._);
          v.m = v.z - midpoint;
        } else {
          v.z = midpoint;
        }
      } else if (w) {
        v.z = w.z + separation(v._, w._);
      }
      v.parent.A = apportion(v, w, v.parent.A || siblings[0]);
    }

    // Computes all real x-coordinates by summing up the modifiers recursively.
    function secondWalk(v) {
      v._.x = v.z + v.parent.m;
      v.m += v.parent.m;
    }

    // The core of the algorithm. Here, a new subtree is combined with the
    // previous subtrees. Threads are used to traverse the inside and outside
    // contours of the left and right subtree up to the highest common level. The
    // vertices used for the traversals are vi+, vi-, vo-, and vo+, where the
    // superscript o means outside and i means inside, the subscript - means left
    // subtree and + means right subtree. For summing up the modifiers along the
    // contour, we use respective variables si+, si-, so-, and so+. Whenever two
    // nodes of the inside contours conflict, we compute the left one of the
    // greatest uncommon ancestors using the function ANCESTOR and call MOVE
    // SUBTREE to shift the subtree and prepare the shifts of smaller subtrees.
    // Finally, we add a new thread (if necessary).
    function apportion(v, w, ancestor) {
      if (w) {
        var vip = v,
            vop = v,
            vim = w,
            vom = vip.parent.children[0],
            sip = vip.m,
            sop = vop.m,
            sim = vim.m,
            som = vom.m,
            shift;
        while (vim = nextRight(vim), vip = nextLeft(vip), vim && vip) {
          vom = nextLeft(vom);
          vop = nextRight(vop);
          vop.a = v;
          shift = vim.z + sim - vip.z - sip + separation(vim._, vip._);
          if (shift > 0) {
            moveSubtree(nextAncestor(vim, v, ancestor), v, shift);
            sip += shift;
            sop += shift;
          }
          sim += vim.m;
          sip += vip.m;
          som += vom.m;
          sop += vop.m;
        }
        if (vim && !nextRight(vop)) {
          vop.t = vim;
          vop.m += sim - sop;
        }
        if (vip && !nextLeft(vom)) {
          vom.t = vip;
          vom.m += sip - som;
          ancestor = v;
        }
      }
      return ancestor;
    }

    function sizeNode(node) {
      node.x *= dx;
      node.y = node.depth * dy;
    }

    tree.separation = function(x) {
      return arguments.length ? (separation = x, tree) : separation;
    };

    tree.size = function(x) {
      return arguments.length ? (nodeSize = false, dx = +x[0], dy = +x[1], tree) : (nodeSize ? null : [dx, dy]);
    };

    tree.nodeSize = function(x) {
      return arguments.length ? (nodeSize = true, dx = +x[0], dy = +x[1], tree) : (nodeSize ? [dx, dy] : null);
    };

    return tree;
  }

  /**
   * This class function creates a TidyTree object.
   * @param {String} newick A valid newick string
   * @param {Object} options A Javascript object containing options to set up the tree
   */
  function TidyTree(data, options, events) {
    let defaults = {
      layout: "vertical",
      type: "tree",
      mode: "smooth",
      leafNodes: true,
      leafLabels: false,
      branchNodes: false,
      branchLabels: false,
      branchDistances: false,
      hStretch: 1,
      vStretch: 1,
      rotation: 0,
      ruler: true,
      animation: 500,
      margin: [50, 50, 50, 50] //CSS order: top, right, bottom, left
    };
    if (!options) options = {};
    Object.assign(this, defaults, options, {
      events: {
        draw: [],
        showtooltip: [],
        hidetooltip: [],
        contextmenu: [],
        search: [],
        select: []
      }
    });

    if(events) Object.keys(events).forEach(e => this.events[e].push(events[e]));

    if (this.parent) this.draw(this.parent);

    if (data instanceof patristic.Branch) {
      this.setData(data);
    } else {
      this.setTree(data);
    }

    if (this.parent) this.recenter();
  }

  /**
   * Update the TidyTree's underlying data structure
   * There are two contexts in which you should call this:
   * 	1. You wish to replace the tree with a completely different tree, given by a different newick string
   * 	2. Your underlying tree data has changed (e.g. the tree has been re-rooted)
   * @param  {Object} data A patristic.Branch object
   * @return {Object}        the TidyTree object
   */
  TidyTree.prototype.setData = function (data) {
    if (!data) throw Error("Invalid Data");
    this.data = data;
    this.range = [Number.MAX_SAFE_INTEGER, Number.MIN_SAFE_INTEGER];
    this.hierarchy = hierarchy(this.data, d => d.children)
      .eachBefore(d => {
        d.value =
          (d.parent ? d.parent.value : 0) + (d.data.length ? d.data.length : 0);
        if (d.value < this.range[0]) this.range[0] = d.value;
        if (d.value > this.range[1]) this.range[1] = d.value;
      })
      .each(d => (d.value /= this.range[1]));
    if (this.parent) return this.redraw();
    return this;
  };

  /**
   * Update the TidyTree's underlying data structure
   * There are two contexts in which you should call this:
   * 	1. You wish to replace the tree with a completely different tree, given by a different newick string
   * 	2. Your underlying tree data has changed (e.g. the tree has been re-rooted)
   * @param  {String} newick A valid newick string
   * @return {Object}        the TidyTree object
   */
  TidyTree.prototype.setTree = function (newick) {
    if (!newick) throw Error("Invalid Newick String");
    return this.setData(patristic.parseNewick(newick));
  };

  /**
   * The available layouts for rendering trees.
   * @type {Array}
   */
  TidyTree.validLayouts = ["horizontal", "vertical", "circular"];

  /**
   * The available types for rendering branches.
   * @type {Array}
   */
  TidyTree.validTypes = ["tree", "weighted", "dendrogram"];

  /**
   * The available modes for rendering branches.
   * @type {Array}
   */
  TidyTree.validModes = ["smooth", "square", "straight"];

  /**
   * Draws a Phylogenetic on the element referred to by selector
   * @param  {String} selector A CSS selector
   * @return {TidyTree}           the TidyTree object
   */
  TidyTree.prototype.draw = function (selector) {
    if (!selector && !this.parent) {
      throw Error("No valid target for drawing given! Where should the tree go?");
    }
    let parent = (this.parent = select(selector ? selector : this.parent));

    this.width =
      parseFloat(parent.style("width")) - this.margin[1] - this.margin[3];
    this.height =
      parseFloat(parent.style("height")) - this.margin[0] - this.margin[2] - 25;

    let tree = d3Tree();

    let svg = parent
      .html(null)
      .append("svg")
      .attr("width", "100%")
      .attr("height", "100%");

    let g = svg.append("g");

    svg
      .append("g")
      .attr("class", "tidytree-ruler")
      .append("rect")
      .attr("y", -5)
      .attr("fill", "white");

    this.zoom = zoom().on("zoom", (event) => {
      let transform = (this.transform = event.transform);
      g.attr(
        "transform",
        `translate(${transform.x},${transform.y}) scale(${transform.k}) rotate(${
        this.rotation
      },${this.layout === "circular" ? 0 : this.width / 2},${
        this.layout === "circular" ? 0 : this.height / 2
      })`
      );
      updateRuler.call(this, transform);
    });
    svg.call(this.zoom);

    g.append("g").attr("class", "tidytree-links");
    g.append("g").attr("class", "tidytree-nodes");

    if (this.events.draw.length) this.events.draw.forEach(c => c());

    return this;
  };

  const getX = d => d.x,
    getY = d => d.y,
    getLength = d => d.weight;

  let linkTransformers = {
    tree: {
      smooth: {
        horizontal: linkHorizontal()
          .x(getY)
          .y(getX),
        vertical: linkVertical()
          .x(getX)
          .y(getY),
        circular: linkRadial()
          .angle(getX)
          .radius(getY)
      },
      straight: {
        horizontal: d =>
          `M${d.source.y} ${d.source.x} L ${d.target.y} ${d.target.x}`,
        vertical: d =>
          `M${d.source.x} ${d.source.y} L ${d.target.x} ${d.target.y}`,
        circular: d => {
          const startAngle = d.source.x - Math.PI / 2,
            startRadius = d.source.y,
            endAngle = d.target.x - Math.PI / 2,
            endRadius = d.target.y;
          const x0 = Math.cos(startAngle),
            y0 = Math.sin(startAngle),
            x1 = Math.cos(endAngle),
            y1 = Math.sin(endAngle);
          return (
            "M" +
            startRadius * x0 +
            "," +
            startRadius * y0 +
            "L" +
            endRadius * x1 +
            "," +
            endRadius * y1
          );
        }
      },
      square: {
        horizontal: d =>
          `M${d.source.y} ${d.source.x} V ${d.target.x} H ${d.target.y}`,
        vertical: d =>
          `M${d.source.x} ${d.source.y} H ${d.target.x} V ${d.target.y}`,
        circular: d => {
          const startAngle = d.source.x - Math.PI / 2,
            startRadius = d.source.y,
            endAngle = d.target.x - Math.PI / 2,
            endRadius = d.target.y;
          const x0 = Math.cos(startAngle),
            y0 = Math.sin(startAngle),
            x1 = Math.cos(endAngle),
            y1 = Math.sin(endAngle);
          return (
            "M" +
            startRadius * x0 +
            "," +
            startRadius * y0 +
            (endAngle === startAngle
              ? ""
              : "A" +
                startRadius +
                "," +
                startRadius +
                " 0 0 " +
                (endAngle > startAngle ? 1 : 0) +
                " " +
                startRadius * x1 +
                "," +
                startRadius * y1) +
            "L" +
            endRadius * x1 +
            "," +
            endRadius * y1
          );
        }
      }
    },
    weighted: {
      smooth: {
        horizontal: linkHorizontal()
          .x(getLength)
          .y(getX),
        vertical: linkVertical()
          .x(getX)
          .y(getLength),
        circular: linkRadial()
          .angle(getX)
          .radius(getLength)
      },
      straight: {
        horizontal: d =>
          `M${d.source.weight} ${d.source.x} L ${d.target.weight} ${d.target.x}`,
        vertical: d =>
          `M${d.source.x} ${d.source.weight} L ${d.target.x} ${d.target.weight}`,
        circular: d => {
          const startAngle = d.source.x - Math.PI / 2,
            startRadius = d.source.weight,
            endAngle = d.target.x - Math.PI / 2,
            endRadius = d.target.weight;
          const x0 = Math.cos(startAngle),
            y0 = Math.sin(startAngle),
            x1 = Math.cos(endAngle),
            y1 = Math.sin(endAngle);
          return (
            "M" +
            startRadius * x0 +
            "," +
            startRadius * y0 +
            "L" +
            endRadius * x1 +
            "," +
            endRadius * y1
          );
        }
      },
      square: {
        horizontal: d =>
          `M${d.source.weight} ${d.source.x} V ${d.target.x} H ${
          d.target.weight
        }`,
        vertical: d =>
          `M${d.source.x} ${d.source.weight} H ${d.target.x} V ${
          d.target.weight
        }`,
        circular: d => {
          const startAngle = d.source.x - Math.PI / 2,
            startRadius = d.source.weight,
            endAngle = d.target.x - Math.PI / 2,
            endRadius = d.target.weight;
          const x0 = Math.cos(startAngle),
            y0 = Math.sin(startAngle),
            x1 = Math.cos(endAngle),
            y1 = Math.sin(endAngle);
          return (
            "M" +
            startRadius * x0 +
            "," +
            startRadius * y0 +
            (endAngle === startAngle
              ? ""
              : "A" +
                startRadius +
                "," +
                startRadius +
                " 0 0 " +
                (endAngle > startAngle ? 1 : 0) +
                " " +
                startRadius * x1 +
                "," +
                startRadius * y1) +
            "L" +
            endRadius * x1 +
            "," +
            endRadius * y1
          );
        }
      }
    }
  };

  linkTransformers.dendrogram = linkTransformers.tree;

  function circularPoint(x, y) {
    return [(y = +y) * Math.cos((x -= Math.PI / 2)), y * Math.sin(x)];
  }

  let nodeTransformers = {
    tree: {
      horizontal: d => `translate(${d.y}, ${d.x})`,
      vertical: d => `translate(${d.x}, ${d.y})`,
      circular: d => `translate(${circularPoint(d.x, d.y)})`
    },
    weighted: {
      horizontal: d => `translate(${d.weight}, ${d.x})`,
      vertical: d => `translate(${d.x}, ${d.weight})`,
      circular: d => `translate(${circularPoint(d.x, d.weight)})`
    }
  };

  nodeTransformers.dendrogram = nodeTransformers.tree;

  const radToDeg = 180 / Math.PI;

  let labelTransformers = {
    tree: {
      straight: {
        horizontal: l =>
          `translate(${(l.source.y + l.target.y) / 2}, ${(l.source.x +
          l.target.x) /
          2}) rotate(${Math.atan(
          (l.target.x - l.source.x) / (l.target.y - l.source.y)
        ) * radToDeg})`,
        vertical: l =>
          `translate(${(l.source.x + l.target.x) / 2}, ${(l.source.y +
          l.target.y) /
          2}) rotate(${Math.atan(
          (l.source.y - l.target.y) / (l.source.x - l.target.x)
        ) * radToDeg})`,
        circular: l => {
          let s = circularPoint(l.source.x, l.source.y),
            t = circularPoint(l.target.x, l.target.y);
          return `translate(${(s[0] + t[0]) / 2}, ${(s[1] + t[1]) /
          2}) rotate(${Math.atan((s[1] - t[1]) / (s[0] - t[0])) * radToDeg})`;
        }
      },
      square: {
        horizontal: l =>
          `translate(${(l.source.y + l.target.y) / 2}, ${l.target.x})`,
        vertical: l =>
          `translate(${l.target.x}, ${(l.source.y + l.target.y) / 2}) rotate(90)`,
        circular: l => {
          let u = circularPoint(l.target.x, (l.source.y + l.target.y) / 2);
          return `translate(${u[0]}, ${u[1]}) rotate(${((l.target.x * radToDeg) %
          180) -
          90})`;
        }
      }
    },
    weighted: {
      straight: {
        horizontal: l =>
          `translate(${(l.source.weight + l.target.weight) / 2}, ${(l.source.x +
          l.target.x) /
          2}) rotate(${Math.atan(
          (l.target.x - l.source.x) / (l.target.weight - l.source.weight)
        ) * radToDeg})`,
        vertical: l =>
          `translate(${(l.source.x + l.target.x) / 2}, ${(l.source.weight +
          l.target.weight) /
          2}) rotate(${Math.atan(
          (l.source.weight - l.target.weight) / (l.source.x - l.target.x)
        ) * radToDeg})`,
        circular: l => {
          let s = circularPoint(l.source.x, l.source.weight),
            t = circularPoint(l.target.x, l.target.weight);
          return `translate(${(s[0] + t[0]) / 2}, ${(s[1] + t[1]) /
          2}) rotate(${Math.atan((s[1] - t[1]) / (s[0] - t[0])) * radToDeg})`;
        }
      },
      square: {
        horizontal: l => `
        translate(${(l.source.weight + l.target.weight) / 2}, ${l.target.x})
      `,
        vertical: l => `
        translate(${l.target.x}, ${(l.source.weight + l.target.weight) / 2})
        rotate(90)
      `,
        circular: l => {
          let u = circularPoint(
            l.target.x,
            (l.source.weight + l.target.weight) / 2
          );
          return `
          translate(${u[0]}, ${u[1]})
          rotate(${((l.target.x * radToDeg) % 180) - 90})
        `;
        }
      }
    }
  };
  labelTransformers.tree.smooth = labelTransformers.tree.straight;
  labelTransformers.weighted.smooth = labelTransformers.weighted.straight;
  labelTransformers.dendrogram = labelTransformers.tree;

  function labeler(d) {
    if (!d.target.data.length) return "0.000";
    return d.target.data.length.toFixed(3);
  }

  /**
   * Redraws the links and relocates the nodes accordingly
   * @return {TidyTree} The TidyTree Object
   */
  TidyTree.prototype.redraw = function () {
    let parent = this.parent;

    this.width  = (parseFloat(parent.style("width" )) - this.margin[1] - this.margin[3]     ) * this.hStretch;
    this.height = (parseFloat(parent.style("height")) - this.margin[0] - this.margin[2] - 25) * this.vStretch;

    this.scalar =
      this.layout === "horizontal" ? this.width :
      this.layout === "vertical" ? this.height :
      Math.min(this.width, this.height) / 2;

    this.hierarchy.each(d => (d.weight = this.scalar * d.value));

    let g = parent.select("svg g");

    let source = (this.type === "tree" ? d3Tree() : cluster()).size(
      this.layout === "circular"   ? [2 * Math.PI, Math.min(this.height, this.width) / 2] :
      this.layout === "horizontal" ? [this.height, this.width] :
      [this.width, this.height]
    );

    if (this.layout === "circular")
      source.separation((a, b) => (a.parent == b.parent ? 1 : 2) / a.depth);

    //Note: You must render links prior to nodes in order to get correct placement!
    let links = g
      .select("g.tidytree-links")
      .selectAll("g.tidytree-link")
      .data(source(this.hierarchy).links(), l => l.source.data._guid + ':' + l.target.data._guid);

    links.join(
      enter => {
        let newLinks = enter.append("g").attr("class", "tidytree-link");

        let linkTransformer = linkTransformers[this.type][this.mode][this.layout];
        newLinks
          .append("path")
          .attr("fill", "none")
          .attr("stroke", "#ccc")
          .attr("d", linkTransformer)
          .transition()
          .duration(this.animation)
          .attr("opacity", 1);

        let labelTransformer = labelTransformers[this.type][this.mode][this.layout];
        newLinks
          .append("text")
          .attr("y", 2)
          .attr("text-anchor", "middle")
          .style("font-size", "12px")
          .text(labeler)
          .attr("transform", labelTransformer)
          .transition()
          .duration(this.animation)
          .style("opacity", this.branchDistances ? 1 : 0);
      },
      update => {
        let linkTransformer = linkTransformers[this.type][this.mode][this.layout];
        let paths = update.select("path");
        if (!this.animation > 0) {
          paths.attr("d", linkTransformer);
        } else {
          paths
            .transition()
            .duration(this.animation / 2)
            .attr("opacity", 0)
            .end()
            .then(() => {
              paths
                .attr("d", linkTransformer)
                .transition()
                .duration(this.animation / 2)
                .attr("opacity", 1);
            });
        }

        let labelTransformer =
          labelTransformers[this.type][this.mode][this.layout];
        let labels = update.select("text");
        if (this.animation) {
          labels
            .transition()
            .duration(this.animation / 2)
            .style("opacity", 0)
            .end()
            .then(() => {
              labels.text(labeler).attr("transform", labelTransformer);
              if (this.branchDistances) {
              labels
                .transition()
                .duration(this.animation / 2)
                .style("opacity", this.branchDistances ? 1 : 0);
              }
            });
        } else {
          labels.text(labeler).attr("transform", labelTransformer);
        }
      },
      exit =>
        exit
          .transition()
          .duration(this.animation)
          .attr("opacity", 0)
          .remove()
    );

    let nodes = g
      .select("g.tidytree-nodes")
      .selectAll("g.tidytree-node")
      .data(this.hierarchy.descendants(), d => d.data._guid);
    nodes.join(
      enter => {
        let nt = nodeTransformers[this.type][this.layout];
        let newNodes = enter
          .append("g")
          .attr("class", "tidytree-node")
          .classed("tidytree-node-internal", d => d.children)
          .classed("tidytree-node-leaf", d => !d.children)
          .attr("transform", nt);

        newNodes
          .append("circle")
          .attr("title", d => d.data.id)
          .style("opacity", d =>
            (d.children && this.branchNodes) ||
            (!d.children && this.leafNodes) ? 1 : 0
          )
          .on("mouseenter focusin", d => this.trigger("showtooltip", d))
          .on("mouseout focusout", d => this.trigger("hidetooltip", d))
          .on("contextmenu", d => this.trigger("contextmenu", d))
          .on("click", d => this.trigger("select", d))
          .attr("r", 2.5);

        let nodeLabels = newNodes
          .append("text")
          .text(d => d.data.id)
          .style("font-size", "12px")
          .attr("y", 2)
          .style("opacity", d =>
            ( d.children && this.branchLabels) ||
            (!d.children && this.leafLabels) ? 1 : 0
          );

        if (this.layout === "vertical") {
          nodeLabels
            .attr("text-anchor", "start")
            .attr("x", 5)
            .transition()
            .duration(this.animation)
            .attr("transform", "rotate(90)");
        } else if (this.layout === "horizontal") {
          nodeLabels
            .attr("text-anchor", "start")
            .attr("x", 5)
            .transition()
            .duration(this.animation)
            .attr("transform", "rotate(0)");
        } else {
          nodeLabels
            .transition()
            .duration(this.animation)
            .attr("transform", l => `rotate(${(((l.x / Math.PI) * 180) % 180) - 90})`)
            .attr("text-anchor", l => l.x % (2 * Math.PI) > Math.PI ? "end" : "start")
            .attr("x", l => (l.x % (2 * Math.PI) > Math.PI ? -5 : 5));
        }

        newNodes
          .transition()
          .duration(this.animation)
          .attr("opacity", 1);
      },
      update => {
        let nodeTransformer = nodeTransformers[this.type][this.layout];
        update
          .transition()
          .duration(this.animation)
          .attr("transform", nodeTransformer);

        let nodeLabels = update.select("text");
        if (this.layout === "vertical") {
          nodeLabels
            .attr("text-anchor", "start")
            .attr("x", 5)
            .transition()
            .duration(this.animation)
            .attr("transform", "rotate(90)");
        } else if (this.layout === "horizontal") {
          nodeLabels
            .attr("text-anchor", "start")
            .attr("x", 5)
            .transition()
            .duration(this.animation)
            .attr("transform", "rotate(0)");
        } else {
          nodeLabels
            .transition()
            .duration(this.animation)
            .attr("transform", l => `rotate(${(((l.x / Math.PI) * 180) % 180) - 90})`)
            .attr("text-anchor", l => l.x % (2 * Math.PI) > Math.PI ? "end" : "start")
            .attr("x", l => (l.x % (2 * Math.PI) > Math.PI ? -5 : 5));
        }
      },
      exit =>
        exit
          .transition()
          .duration(this.animation)
          .attr("opacity", 0)
          .remove()
    );

    updateRuler.call(this);

    return this;
  };

  function updateRuler(transform) {
    if (!transform) transform = { k: 1 };
    let height = parseFloat(this.parent.style("height")) - this.margin[2] - 15;
    let ruler = this.parent.select("g.tidytree-ruler");
    let bg = ruler.select("rect");
    if (this.ruler) {
      if (this.layout == "horizontal") {
        ruler.attr("transform", `translate(${this.margin[3]}, ${height})`);
        bg
          .attr("width", `calc(100% - ${this.margin[1] + this.margin[3] - 15}px)`)
          .attr("height", "25px")
          .attr("x", -5);
      } else {
        ruler.attr("transform", `translate(${this.margin[3] - 10}, ${this.margin[0]})`);
        bg
          .attr("height", `calc(100% - ${this.margin[0] + this.margin[2] - 15}px)`)
          .attr("width", "25px")
          .attr("x", -25);
      }
      let axis = this.layout == "horizontal" ? axisBottom() : axisLeft();
      if (this.type === "tree" && this.layout !== "circular") {
        ruler
          .attr("opacity", 1)
          .call(
            axis.scale(
              linear$1(
                [0, this.hierarchy.height / transform.k],
                [0, this.scalar]
              )
            )
          );
      } else if (this.type === "weighted" && this.layout !== "circular") {
        ruler
          .attr("opacity", 1)
          .call(
            axis.scale(
              linear$1(
                [this.range[0], this.range[1] / transform.k],
                [0, this.scalar]
              )
            )
          );
      } else {
        ruler
          .transition()
          .duration(this.animation)
          .attr("opacity", 0);
      }
    } else {
      ruler
        .transition()
        .duration(this.animation)
        .attr("opacity", 0);
    }
  }

  /**
   * Recenters the tree in the center of the view
   * @return {TidyTree} The TidyTree object
   */
  TidyTree.prototype.recenter = function () {
    let svg = this.parent.select("svg"),
      x = this.margin[0],
      y = this.margin[3];
    if (this.layout === "circular") {
      x += parseFloat(svg.style("width")) / 2;
      y += parseFloat(svg.style("height")) / 2;
    }
    svg
      .transition()
      .duration(this.animation)
      .call(this.zoom.transform, identity$4.translate(x, y));
    return this;
  };

  /**
   * Set the TidyTree's layout
   * @param {String} newLayout The new layout
   * @return {TidyTree} The TidyTree Object
   */
  TidyTree.prototype.setLayout = function (newLayout) {
    if (!TidyTree.validLayouts.includes(newLayout)) {
      throw Error(`
      Cannot set TidyTree to layout: ${newLayout}\n
      Valid layouts are: ${TidyTree.validLayouts.join(', ')}
    `);
    }
    this.layout = newLayout;
    if (this.parent) return this.redraw();
    return this;
  };

  /**
   * Set the TidyTree's mode
   * @param {String} newMode The new mode
   * @return {TidyTree} The TidyTree object
   */
  TidyTree.prototype.setMode = function (newMode) {
    if (!TidyTree.validModes.includes(newMode)) {
      throw Error(`
      Cannot set TidyTree to mode: ${newMode},\n
      Valid modes are: ${TidyTree.validModes.join(', ')}
    `);
    }
    this.mode = newMode;
    if (this.parent) return this.redraw();
    return this;
  };

  /**
   * Set the TidyTree's type
   * @param {Boolean} newType The new type
   * @return {TidyTree} the TidyTree object
   */
  TidyTree.prototype.setType = function (newType) {
    if (!TidyTree.validTypes.includes(newType)) {
      throw Error(`
      Cannot set TidyTree to type: ${newType},\n
      Valid types are: ${TidyTree.validTypes.join(', ')}
    `);
    }
    this.type = newType;
    if (this.parent) return this.redraw();
    return this;
  };

  /**
   * Set the TidyTree's rotation
   * @param {Number} degrees The new number of degrees by which to rotate the tree
   * @return {TidyTree} the TidyTree object
   */
  TidyTree.prototype.setRotation = function (degrees) {
    this.rotation = degrees;
    if (this.parent)
      this.parent
        .select("svg g")
        .attr("transform", `
        translate(${this.transform.x},${this.transform.y})
        scale(${this.transform.k})
        rotate(${this.rotation},
          ${this.layout === "circular" ? 0 : this.width / 2},
          ${this.layout === "circular" ? 0 : this.height / 2}
        )
      `);
    return this;
  };

  /**
   * Set the TidyTree's Horizontal Stretch
   * @param {Number} proportion The new proportion by which to stretch the tree
   * @return {TidyTree} the TidyTree object
   */
  TidyTree.prototype.setHStretch = function (proportion) {
    this.hStretch = parseFloat(proportion);
    if (this.parent) {
      let animCache = this.animation;
      this.setAnimation(0);
      this.redraw();
      this.setAnimation(animCache);
    }
    return this;
  };

  /**
   * Set the TidyTree's Vertical Stretch
   * @param {Number} proportion The new proportion by which to stretch the tree
   * @return {TidyTree} the TidyTree object
   */
  TidyTree.prototype.setVStretch = function (proportion) {
    this.vStretch = parseFloat(proportion);
    if (this.parent) {
      let animCache = this.animation;
      this.setAnimation(0);
      this.redraw();
      this.setAnimation(animCache);
    }
    return this;
  };

  /**
   * Set the TidyTree's animation time. Note that this does not trigger a
   * redraw.
   * @param {number} time The desired duration of an animation, in ms. Set to 0
   * to turn animations off completely.
   * @return {TidyTree} The TidyTree object
   */
  TidyTree.prototype.setAnimation = function (time) {
    this.animation = time;
    return this;
  };

  /**
   * Shows or hides the Branch Nodes
   * @param  {Boolean} show Should Branch nodes be shown?
   * @return {TidyTree} the TidyTree object
   */
  TidyTree.prototype.setBranchNodes = function (show) {
    this.branchNodes = show ? true : false;
    if (this.parent) {
      //i.e. has already been drawn
      this.parent
        .select("svg")
        .selectAll("g.tidytree-node-internal circle")
        .transition()
        .duration(this.animation)
        .style("opacity", show ? 1 : 0);
    }
    return this;
  };

  /**
   * Restyles Leaf Nodes
   * @param  {Function} styler A function that restyles each node. `styler`
   * receives a reference to the DOM node to be styled, and an associated data
   * object.
   * @return {TidyTree} the TidyTree Object
   */
  TidyTree.prototype.eachBranchNode = function (styler) {
    if (!this.parent)
      throw Error(
        "Tree has not been rendered yet! Can't style Nodes that don't exist!"
      );
    this.parent
      .select("svg")
      .selectAll("g.tidytree-node-internal circle")
      .each(function (d) { styler(this, d); });
    return this;
  };

  /**
   * Set the TidyTree's branchLabels
   * @param  {Boolean} show Should the TidyTree show branchLabels?
   * @return {TidyTree}     the TidyTree Object
   */
  TidyTree.prototype.setBranchLabels = function (show) {
    this.branchLabels = show ? true : false;
    if (this.parent) {
      //i.e. has already been drawn
      this.parent
        .select("svg")
        .selectAll("g.tidytree-node-internal text")
        .transition()
        .duration(this.animation)
        .style("opacity", show ? 1 : 0);
    }
    return this;
  };

  /**
   * Restyles Branch Label
   * @param  {Function} styler A function that restyles each node. `styler`
   * receives a reference to the DOM node to be styled, and an associated data
   * object.
   * @return {TidyTree} the TidyTree Object
   */
  TidyTree.prototype.eachBranchLabel = function (styler) {
    if (!this.parent){
      throw Error("Tree has not been rendered yet! Can't style Nodes that don't exist!");
    }
    this.parent
      .select("svg")
      .selectAll("g.tidytree-node-internal text")
      .each(function (d, i, l) { styler(this, d); });
    return this;
  };

  /**
   * Shows or hides the TidyTree's branch labels
   * @param {Boolean} show Should the TidyTree show branch distances?
   * @return {TidyTree} The TidyTree Object
   */
  TidyTree.prototype.setBranchDistances = function (show) {
    this.branchDistances = show ? true : false;
    if (this.parent) {
      //i.e. has already been drawn
      let links = this.parent
        .select("svg g.tidytree-links")
        .selectAll("g.tidytree-link")
        .selectAll("text");
      links.attr("transform", labelTransformers[this.type][this.mode][this.layout]);
      links
        .transition()
        .duration(this.animation)
        .style("opacity", show ? 1 : 0);
    }
    return this;
  };

  /**
   * Restyles Branch Distances
   * @param  {Function} styler A function that restyles each node. `styler`
   * receives a reference to the DOM node to be styled, and an associated data
   * object.
   * @return {TidyTree} the TidyTree Object
   */
  TidyTree.prototype.eachBranchDistance = function (styler) {
    if (!this.parent)
      throw Error("Tree has not been rendered yet! Can't style Nodes that don't exist!");
    this.parent
      .select("svg g.tidytree-links")
      .selectAll("g.tidytree-link")
      .selectAll("text")
      .each(function (d, i, l) { styler(this, d); });
    return this;
  };

  /**
   * Shows or Hides the Leaf Nodes
   * @param  {Boolean} show Should leaf nodes be visible?
   * @return {TidyTree} The TidyTree Object
   */
  TidyTree.prototype.setLeafNodes = function (show) {
    this.leafNodes = show ? true : false;
    if (this.parent) {
      //i.e. has already been drawn
      this.parent
        .select("svg")
        .selectAll("g.tidytree-node-leaf circle")
        .transition()
        .duration(this.animation)
        .style("opacity", show ? 1 : 0);
    }
    return this;
  };

  /**
   * Restyles Leaf Nodes
   * @param  {Function} styler A function that restyles each node. `styler`
   * receives a reference to the DOM node to be styled, and an associated data
   * object.
   * @return {TidyTree} the TidyTree Object
   */
  TidyTree.prototype.eachLeafNode = function (styler) {
    if (!this.parent){
      throw Error("Tree has not been rendered yet! Can't style Nodes that don't exist!");
    }
    this.parent
      .select("svg")
      .selectAll("g.tidytree-node-leaf circle")
      .each(function (d) {
        styler(this, d);
      });
    return this;
  };

  /**
   * Shows or Hides the TidyTree's Leaf Labels
   * @param  {Boolean} show Should the TidyTree show leafLabels?
   * @return {TidyTree}     the TidyTree Object
   */
  TidyTree.prototype.setLeafLabels = function (show) {
    this.leafLabels = show ? true : false;
    if (this.parent) {
      //i.e. has already been drawn
      this.parent
        .select("svg")
        .selectAll("g.tidytree-node-leaf text")
        .transition()
        .duration(this.animation)
        .style("opacity", show ? 1 : 0);
    }
    return this;
  };

  /**
   * Restyles Leaf Labels
   * @param  {Function} styler A function that restyles each node. `styler`
   * receives a reference to the DOM node to be styled, and an associated data
   * object.
   * @return {TidyTree} the TidyTree Object
   */
  TidyTree.prototype.eachLeafLabel = function (styler) {
    if (!this.parent){
      throw Error("Tree has not been rendered yet! Can't style Nodes that don't exist!");
    }
    this.parent
      .select("svg")
      .selectAll("g.tidytree-node-leaf text")
      .each(function (d) { styler(this, d); });
    return this;
  };

  /**
   * Shows or hides the TidyTree's branch labels
   * @param {Boolean} show Should the TidyTree show branchLabels?
   * @return {TidyTree} The TidyTree Object
   */
  TidyTree.prototype.setRuler = function (show) {
    this.ruler = show ? true : false;
    if (this.parent) {
      //i.e. has already been drawn
      if (show) {
        this.parent
          .select("g.tidytree-ruler")
          .transition()
          .duration(this.animation)
          .attr("opacity", 1);
      } else {
        this.parent
          .select("g.tidytree-ruler")
          .transition()
          .duration(this.animation)
          .attr("opacity", 0);
      }
    }
    return this;
  };

  /**
   * Searches the tree, returns Search Results
   * @param  {Function} test A function which takes a Branch and returns a Truthy
   * or Falsy value.
   * @return {Array} The array of results
   */
  TidyTree.prototype.search = function (test) {
    if (!test) return;
    let results = this.parent
      .select("svg g.tidytree-nodes")
      .selectAll("g.tidytree-node")
      .filter(test);
    if (this.events.search.length) this.events.search.forEach(c => c(results));
    return results;
  };

  /**
   * Attaches a new event listener
   * Please note that this is not yet functioning.
   * @param  {String}   events   A space-delimited list of event names
   * @param  {Function} callback The function to run when one of the `events` occurs.
   * @return {TidyTree} The TidyTree on which this method was called.
   */
  TidyTree.prototype.on = function (events, callback) {
    events.split(" ").forEach(event => this.events[event].push(callback));
    return this;
  };

  /**
   * Removes all event listeners from the given events
   * @param  {String}   events   A space-delimited list of event names
   * @return {TidyTree} The TidyTree on which this method was called.
   */
  TidyTree.prototype.off = function (events) {
    events.split(" ").forEach(event => (this.events[event] = []));
    return this;
  };

  /**
   * Forces the tree to respond as though an `event` has occurred
   * @param  {String} events space-delimited list of names of events to trigger.
   * @param  {Spread} args Any arguments which should be passed to the event
   * handler(s).
   * @return The output of the callback run on `event`
   */
  TidyTree.prototype.trigger = function (events, ...args) {
    return events.split(" ").map(event => {
      if (this.events[event].length)
        return this.events[event].map(handler => handler(args));
      return [];
    });
  };

  /**
   * Destroys the TidyTree
   * @return {undefined}
   */
  TidyTree.prototype.destroy = function () {
    if (this.parent) {
      //i.e. has already been drawn
      this.parent.html(null);
    }
  };

  return TidyTree;

}());
